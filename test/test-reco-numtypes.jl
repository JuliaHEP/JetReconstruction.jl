# Test reconstruction with different numerical types

include("common.jl")

for TestType in [Float32, Float16]
    let inputs_ee = JetReconstruction.read_final_state_particles(events_file_ee,
                                                                 EEJet{TestType};
                                                                 maxevents = 1)
        @testset "EE jet reconstruction with $TestType" begin
            # EE Algorithms
            @test typeof(inputs_ee[1]) == Vector{EEJet{TestType}}
            @test typeof(jet_reconstruct(inputs_ee[1]; algorithm = JetAlgorithm.Durham)) ==
                  ClusterSequence{TestType, EEJet{TestType}}
            @test typeof(jet_reconstruct(inputs_ee[1]; algorithm = JetAlgorithm.EEKt,
                                         p = -1,
                                         R = 1.0)) ==
                  ClusterSequence{TestType, EEJet{TestType}}
        end
    end
    let inputs_pp = JetReconstruction.read_final_state_particles(events_file_pp,
                                                                 PseudoJet{TestType};
                                                                 maxevents = 1)
        @testset "PP jet reconstruction with $TestType" begin
            # PP Algorithms
            @test typeof(inputs_pp[1]) == Vector{PseudoJet{TestType}}
            @test typeof(jet_reconstruct(inputs_pp[1]; algorithm = JetAlgorithm.AntiKt)) ==
                  ClusterSequence{TestType, PseudoJet{TestType}}
            @test typeof(jet_reconstruct(inputs_pp[1]; algorithm = JetAlgorithm.CA)) ==
                  ClusterSequence{TestType, PseudoJet{TestType}}
        end
    end
end
