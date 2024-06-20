#! /usr/bin/env julia

using JetReconstruction
using Test
using JSON
using LorentzVectorHEP
using Logging

const events_file = joinpath(@__DIR__, "data", "events.hepmc3.gz")

const algorithms = Dict(-1 => "Anti-kt",
                        0 => "Cambridge/Achen",
                        1 => "Inclusive-kt")

"""Simple structure with necessary parameters for an exclusive selection test"""
struct InclusiveTest
    algname::AbstractString
    selction::AbstractString
    power::Int
    fastjet_file::AbstractString
    dijmax::Any
    njets::Any
    function InclusiveTest(selction, power, fastjet_file, dijmax, njets)
        new(algorithms[power], selction, power, fastjet_file, dijmax, njets)
    end
end

"""Read JSON file with fastjet jets in it"""
function read_fastjet_outputs(fname)
    f = open(fname)
    JSON.parse(f)
end

"""Sort jet outputs by pt of final jets"""
function sort_jets!(event_jet_array)
    jet_pt(jet) = jet["pt"]
    for e in event_jet_array
        sort!(e["jets"], by = jet_pt, rev = true)
    end
end

function sort_jets!(jet_array::Vector{FinalJet})
    jet_pt(jet) = jet.pt
    sort!(jet_array, by = jet_pt, rev = true)
end

function sort_jets!(jet_array::Vector{LorentzVectorCyl})
    jet_pt(jet) = jet.pt
    sort!(jet_array, by = jet_pt, rev = true)
end

function main()
    # Read our fastjet inclusive outputs (we read for anti-kt, cambridge/achen, inclusive-kt)
    fastjet_alg_files_inclusive = Dict(-1 => joinpath(@__DIR__, "data",
                                                      "jet-collections-fastjet-inclusive-akt.json"),
                                       0 => joinpath(@__DIR__, "data",
                                                     "jet-collections-fastjet-inclusive-ca.json"),
                                       1 => joinpath(@__DIR__, "data",
                                                     "jet-collections-fastjet-inclusive-ikt.json"))

    fastjet_data = Dict{Int, Vector{Any}}()
    for (k, v) in pairs(fastjet_alg_files_inclusive)
        fastjet_jets = read_fastjet_outputs(v)
        sort_jets!(fastjet_jets)
        fastjet_data[k] = fastjet_jets
    end

    # Test each stratgy for inclusive jet selection
    for power in keys(algorithms)
        do_test_compare_to_fastjet(RecoStrategy.N2Plain, fastjet_data[power],
                                   algname = algorithms[power], power = power)
        do_test_compare_to_fastjet(RecoStrategy.N2Tiled, fastjet_data[power],
                                   algname = algorithms[power], power = power)
    end

    # Compare inputing data in PseudoJet with using a LorentzVector
    do_test_compare_types(RecoStrategy.N2Plain, algname = algorithms[-1], power = -1)
    do_test_compare_types(RecoStrategy.N2Tiled, algname = algorithms[-1], power = -1)

    # Now test exclusive selections
    inclusive_tests = InclusiveTest[]
    push!(inclusive_tests,
          InclusiveTest("exclusive njets=4", 1, "jet-collections-fastjet-njets4-ikt.json",
                        nothing, 4))
    push!(inclusive_tests,
          InclusiveTest("exclusive dijmax=20", 1, "jet-collections-fastjet-dij20-ikt.json",
                        20.0, nothing))
    push!(inclusive_tests,
          InclusiveTest("exclusive njets=4", 0, "jet-collections-fastjet-njets4-ca.json",
                        nothing, 4))
    push!(inclusive_tests,
          InclusiveTest("exclusive dijmax=0.99", 0,
                        "jet-collections-fastjet-dij099-ca.json", 0.99, nothing))

    for test in inclusive_tests
        input_file = joinpath(@__DIR__, "data", test.fastjet_file)
        fastjet_jets = read_fastjet_outputs(input_file)
        sort_jets!(fastjet_jets)
        do_test_compare_to_fastjet(RecoStrategy.N2Tiled, fastjet_jets;
                                   algname = test.algname,
                                   power = test.power, selection = test.selction,
                                   njets = test.njets, dijmax = test.dijmax)
    end
end

"""
Run the jet test

dijmax -> apply inclusive dijmax cut
njets -> apply inclusive njets cut

If dijmax and njets are nothing, test inclusive jets with pt >= ptmin
"""
function do_test_compare_to_fastjet(strategy::RecoStrategy.Strategy, fastjet_jets;
                                    algname = "Unknown",
                                    selection = "Inclusive",
                                    ptmin::Float64 = 5.0,
                                    distance::Float64 = 0.4,
                                    power::Integer = -1,
                                    dijmax = nothing,
                                    njets = nothing)

    # Strategy
    if (strategy == RecoStrategy.N2Plain)
        jet_reconstruction = plain_jet_reconstruct
        strategy_name = "N2Plain"
    elseif (strategy == RecoStrategy.N2Tiled)
        jet_reconstruction = tiled_jet_reconstruct
        strategy_name = "N2Tiled"
    elseif (strategy == RecoStrategy.Best)
        jet_reconstruction = jet_reconstruct
        strategy_name = "Best"
    else
        throw(ErrorException("Strategy not yet implemented"))
    end

    # Now run our jet reconstruction...
    # From PseudoJets
    events::Vector{Vector{PseudoJet}} = read_final_state_particles(events_file)
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        cluster_seq = jet_reconstruction(event, R = distance, p = power)
        if !isnothing(dijmax)
            finaljets = final_jets(exclusive_jets(cluster_seq, dcut = dijmax))
        elseif !isnothing(njets)
            finaljets = final_jets(exclusive_jets(cluster_seq, njets = njets))
        else
            finaljets = final_jets(inclusive_jets(cluster_seq, ptmin))
        end
        sort_jets!(finaljets)
        push!(jet_collection, FinalJets(ievt, finaljets))
    end

    @testset "Jet Reconstruction compare to FastJet: Strategy $strategy_name, Algorithm $algname, Selection $selection" begin
        # Test each event in turn...
        for (ievt, event) in enumerate(jet_collection)
            @testset "Event $(ievt)" begin
                @test size(event.jets) == size(fastjet_jets[ievt]["jets"])
                # Test each jet in turn
                for (ijet, jet) in enumerate(event.jets)
                    if ijet <= size(fastjet_jets[ievt]["jets"])[1]
                        # Approximate test - note that @test macro passes the 
                        # tolerance into the isapprox() function
                        # Use atol for position as this is absolute, but rtol for
                        # the momentum
                        # Sometimes phi could be in the range [-π, π], but FastJet always is [0, 2π]
                        normalised_phi = jet.phi < 0.0 ? jet.phi + 2π : jet.phi
                        @test jet.rap≈fastjet_jets[ievt]["jets"][ijet]["rap"] atol=1e-7
                        @test normalised_phi≈fastjet_jets[ievt]["jets"][ijet]["phi"] atol=1e-7
                        @test jet.pt≈fastjet_jets[ievt]["jets"][ijet]["pt"] rtol=1e-6
                    end
                end
            end
        end
    end
end

function do_test_compare_types(strategy::RecoStrategy.Strategy;
                               algname = "Unknown",
                               ptmin::Float64 = 5.0,
                               distance::Float64 = 0.4,
                               power::Integer = -1)

    # Strategy
    if (strategy == RecoStrategy.N2Plain)
        jet_reconstruction = plain_jet_reconstruct
        strategy_name = "N2Plain"
    elseif (strategy == RecoStrategy.N2Tiled)
        jet_reconstruction = tiled_jet_reconstruct
        strategy_name = "N2Tiled"
    else
        throw(ErrorException("Strategy not yet implemented"))
    end

    # Now run our jet reconstruction...
    # From PseudoJets
    events::Vector{Vector{PseudoJet}} = read_final_state_particles(events_file)
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        finaljets = final_jets(inclusive_jets(jet_reconstruction(event, R = distance,
                                                                 p = power), ptmin))
        sort_jets!(finaljets)
        push!(jet_collection, FinalJets(ievt, finaljets))
    end

    # From LorentzVector
    events_lv::Vector{Vector{LorentzVector}} = read_final_state_particles_lv(events_file)
    jet_collection_lv = FinalJets[]
    for (ievt, event) in enumerate(events_lv)
        finaljets = final_jets(inclusive_jets(jet_reconstruction(event, R = distance,
                                                                 p = power), ptmin))
        sort_jets!(finaljets)
        push!(jet_collection_lv, FinalJets(ievt, finaljets))
    end

    @testset "Jet Reconstruction Compare PseudoJet and LorentzVector, Strategy $strategy_name, Algorithm $algname" begin
        # Here we test that inputing LorentzVector gave the same results as PseudoJets
        for (ievt, (event, event_lv)) in enumerate(zip(jet_collection, jet_collection_lv))
            @testset "Event $(ievt)" begin
                @test size(event.jets) == size(event_lv.jets)
                # Test each jet in turn
                for (jet, jet_lv) in zip(event.jets, event_lv.jets)
                    @test jet.rap≈jet_lv.rap atol=1e-7
                    @test jet.phi≈jet_lv.phi atol=1e-7
                    @test jet.pt≈jet_lv.pt rtol=1e-6
                end
            end
        end
    end
end

logger = ConsoleLogger(stdout, Logging.Warn)
global_logger(logger)
main()
