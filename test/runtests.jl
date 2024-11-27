#! /usr/bin/env julia

using JetReconstruction
using Test
using JSON
using LorentzVectorHEP
using Logging

include("common.jl")

function main()
    # A few unit tests
    @testset "Algorithm/power consistency" begin
        @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.AntiKt,
                                                                  p = -1)
        @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.CA,
                                                                  p = 0)
        @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.Kt,
                                                                  p = 1)

        @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.AntiKt,
                                                                  p = nothing)
        @test JetReconstruction.check_algorithm_power_consistency(algorithm = nothing,
                                                                  p = -1)

        @test_throws ArgumentError JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.AntiKt,
                                                                                       p = 0)
        @test_throws ArgumentError JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.Kt,
                                                                                       p = 1.5)

        @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.GenKt,
                                                                  p = 1.5)
        @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.GenKt,
                                                                  p = -0.5)

        @test_throws ArgumentError JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.GenKt,
                                                                                       p = nothing)
    end

    # New test structure, factorised tests for pp and e+e-
    include("test-pp-reconstruction.jl")
    include("test-ee-reconstruction.jl")

    # Compare inputting data in PseudoJet with using a LorentzVector
    do_test_compare_types(RecoStrategy.N2Plain, algname = pp_algorithms[-1], power = -1)
    do_test_compare_types(RecoStrategy.N2Tiled, algname = pp_algorithms[-1], power = -1)

    # Check jet constituents
    include("test-constituents.jl")

    # Suppress these tests for now, as the examples Project.toml is rather heavy
    # because of the GLMakie dependency, plus on a CI there is no GL subsystem,
    # so things fail. The examples should be restructured to have a cleaner set
    # of examples that are good to run in the CI.
    #
    # Now run a few tests with our examples
    # include("tests_examples.jl")
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
                                    power::Real = -1,
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
    events::Vector{Vector{PseudoJet}} = read_final_state_particles(events_file_pp)
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        # First run the reconstruction
        if power == 1.5
            cluster_seq = jet_reconstruction(event, R = distance,
                                             algorithm = JetAlgorithm.GenKt, p = power)
        else
            cluster_seq = jet_reconstruction(event, R = distance, p = power)
        end
        # Now make the requested selection
        if !isnothing(dijmax)
            selected_jets = exclusive_jets(cluster_seq; dcut = dijmax)
        elseif !isnothing(njets)
            selected_jets = exclusive_jets(cluster_seq; njets = njets)
        else
            selected_jets = inclusive_jets(cluster_seq; ptmin = ptmin)
        end
        # And extract in out final_jets format
        finaljets = final_jets(selected_jets)
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
    events::Vector{Vector{PseudoJet}} = read_final_state_particles(events_file_pp)
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        finaljets = final_jets(inclusive_jets(jet_reconstruction(event, R = distance,
                                                                 p = power), ptmin = ptmin))
        sort_jets!(finaljets)
        push!(jet_collection, FinalJets(ievt, finaljets))
    end

    # From LorentzVector
    events_lv::Vector{Vector{LorentzVector}} = read_final_state_particles(events_file_pp;
                                                                          T = LorentzVector)
    jet_collection_lv = FinalJets[]
    for (ievt, event) in enumerate(events_lv)
        finaljets = final_jets(inclusive_jets(jet_reconstruction(event, R = distance,
                                                                 p = power), ptmin = ptmin))
        sort_jets!(finaljets)
        push!(jet_collection_lv, FinalJets(ievt, finaljets))
    end

    @testset "Jet Reconstruction Compare PseudoJet and LorentzVector, Strategy $strategy_name, Algorithm $algname" begin
        # Here we test that inputting LorentzVector gave the same results as PseudoJets
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
