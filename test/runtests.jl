#! /usr/bin/env julia

using JetReconstruction
using Test
using JSON
using LorentzVectorHEP
using Logging

"""Read JSON file with fastjet jets in it"""
function read_fastjet_outputs(; fname = "test/data/jet_collections_fastjet_akt.json")
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

function main()
    # Read our fastjet outputs (we read for anti-kt, cambridge/achen, inclusive-kt)
    algorithms = Dict(-1 => "Anti-kt",
        0 => "Cambridge/Achen",
        1 => "Inclusive-kt")
    fastjet_alg_files = Dict(-1 => "test/data/jet_collections_fastjet_akt.json",
        0 => "test/data/jet_collections_fastjet_ca.json",
        1 => "test/data/jet_collections_fastjet_ikt.json")

    fastjet_data = Dict{Int, Vector{Any}}()
    for (k, v) in pairs(fastjet_alg_files)
        fastjet_jets = read_fastjet_outputs(fname = v)
        sort_jets!(fastjet_jets)
        fastjet_data[k] = fastjet_jets
    end

    # Test each stratgy...
    for power in keys(algorithms)
        do_test_compare_to_fastjet(N2Plain, fastjet_data[power], algname = algorithms[power], power = power)
        do_test_compare_to_fastjet(N2Tiled, fastjet_data[power], algname = algorithms[power], power = power)
    end

    # Compare inputing data in PseudoJet with using a LorentzVector
    do_test_compare_types(N2Plain, algname = algorithms[-1], power = -1)
    do_test_compare_types(N2Tiled, algname = algorithms[-1], power = -1)
end

function do_test_compare_to_fastjet(strategy::JetRecoStrategy, fastjet_jets;
    algname = "Unknown",
    ptmin::Float64 = 5.0,
    distance::Float64 = 0.4,
    power::Integer = -1)

    # Strategy
    if (strategy == N2Plain)
        jet_reconstruction = plain_jet_reconstruct
        strategy_name = "N2Plain"
    elseif (strategy == N2Tiled)
        jet_reconstruction = tiled_jet_reconstruct
        strategy_name = "N2Tiled"
    else
        throw(ErrorException("Strategy not yet implemented"))
    end

    # Now run our jet reconstruction...
    # From PseudoJets
    events::Vector{Vector{PseudoJet}} = read_final_state_particles("test/data/events.hepmc3")
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        finaljets, _ = jet_reconstruction(event, R = distance, p = power, ptmin = ptmin)
        fj = final_jets(finaljets, ptmin)
        sort_jets!(fj)
        push!(jet_collection, FinalJets(ievt, fj))
    end

    @testset "Jet Reconstruction compare to FastJet: Strategy $strategy_name, Algorithm $algname" begin
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
                        @test jet.rap ≈ fastjet_jets[ievt]["jets"][ijet]["rap"] atol = 1e-7
                        @test normalised_phi ≈ fastjet_jets[ievt]["jets"][ijet]["phi"] atol = 1e-7
                        @test jet.pt ≈ fastjet_jets[ievt]["jets"][ijet]["pt"] rtol = 1e-6
                    end
                end
            end
        end
    end
end

function do_test_compare_types(strategy::JetRecoStrategy;
    algname = "Unknown",
    ptmin::Float64 = 5.0,
    distance::Float64 = 0.4,
    power::Integer = -1)

    # Strategy
    if (strategy == N2Plain)
        jet_reconstruction = plain_jet_reconstruct
        strategy_name = "N2Plain"
    elseif (strategy == N2Tiled)
        jet_reconstruction = tiled_jet_reconstruct
        strategy_name = "N2Tiled"
    else
        throw(ErrorException("Strategy not yet implemented"))
    end

    # Now run our jet reconstruction...
    # From PseudoJets
    events::Vector{Vector{PseudoJet}} = read_final_state_particles("test/data/events.hepmc3")
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        finaljets, _ = jet_reconstruction(event, R = distance, p = power, ptmin = ptmin)
        fj = final_jets(finaljets, ptmin)
        sort_jets!(fj)
        push!(jet_collection, FinalJets(ievt, fj))
    end

    # From LorentzVector
    events_lv::Vector{Vector{LorentzVector}} = read_final_state_particles_lv("test/data/events.hepmc3")
    jet_collection_lv = FinalJets[]
    for (ievt, event) in enumerate(events_lv)
        finaljets, _ = jet_reconstruction(event, R = distance, p = power, ptmin = ptmin)
        fj = final_jets(finaljets, ptmin)
        sort_jets!(fj)
        push!(jet_collection_lv, FinalJets(ievt, fj))
    end

    @testset "Jet Reconstruction Compare PseudoJet and LorentzVector, Strategy $strategy_name, Algorithm $algname" begin
        # Here we test that inputing LorentzVector gave the same results as PseudoJets
        for (ievt, (event, event_lv)) in enumerate(zip(jet_collection, jet_collection_lv))
            @testset "Event $(ievt)" begin
                @test size(event.jets) == size(event_lv.jets)
                # Test each jet in turn
                for (jet, jet_lv) in zip(event.jets, event_lv.jets)
                    @test jet.rap ≈ jet_lv.rap atol = 1e-7
                    @test jet.phi ≈ jet_lv.phi atol = 1e-7
                    @test jet.pt ≈ jet_lv.pt rtol = 1e-6
                end
            end
        end
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    logger = ConsoleLogger(stdout, Logging.Warn)
    global_logger(logger)
    main()
end
