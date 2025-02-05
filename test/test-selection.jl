# Tests of jets selection functionality

include("common.jl")

""" Helper function to test exclusive jets selection with given algorithm power argument."""
function test_exclusive_jets_power(event, algo; p = nothing, should_throw::Bool)
    test_cone_size = 0.4
    test_strategy = RecoStrategy.Best
    test_njets = 4

    @testset "Exclusive jets selection: alg=$(algo) p=$(p)" begin
        cluster_seq = JetReconstruction.jet_reconstruct(event; algorithm = algo, p = p,
                                                        R = test_cone_size,
                                                        strategy = test_strategy)
        if should_throw
            @test_throws ArgumentError begin
                JetReconstruction.exclusive_jets(cluster_seq; njets = test_njets)
            end
        else
            @test JetReconstruction.exclusive_jets(cluster_seq; njets = test_njets) isa Any
        end
    end
end

@testset "Exclusive jets selection algorithm power requirements" begin
    pp_event = first(JetReconstruction.read_final_state_particles(events_file_pp))
    ee_event = first(JetReconstruction.read_final_state_particles(events_file_ee))

    # Tests for pp_event
    test_exclusive_jets_power(pp_event, JetAlgorithm.GenKt; p = -1.0,
                              should_throw = true)
    test_exclusive_jets_power(pp_event, JetAlgorithm.GenKt; p = 1.0,
                              should_throw = false)
    test_exclusive_jets_power(pp_event, JetAlgorithm.AntiKt; should_throw = true)
    for alg in [JetAlgorithm.CA, JetAlgorithm.Kt]
        test_exclusive_jets_power(pp_event, alg; should_throw = false)
    end

    # Tests for ee_event
    test_exclusive_jets_power(ee_event, JetAlgorithm.EEKt; p = -1.0,
                              should_throw = true)
    test_exclusive_jets_power(ee_event, JetAlgorithm.EEKt; p = 1.0,
                              should_throw = false)
    test_exclusive_jets_power(ee_event, JetAlgorithm.Durham; should_throw = false)
end
