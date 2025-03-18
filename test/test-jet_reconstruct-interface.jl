# Test the different interface modes of the main jet_reconstruct() function

include("common.jl")

let inputs = JetReconstruction.read_final_state_particles(events_file_ee)
    @testset "jet_reconstruct() interface" begin
        # EE Algorithms
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.Durham)) ==
              ClusterSequence{EEJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.EEKt, p = -1,
                                     R = 1.0)) == ClusterSequence{EEJet}
        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.EEKt)

        # PP Algorithms
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.AntiKt)) ==
              ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.CA)) ==
              ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.Kt)) ==
              ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.AntiKt, p = -1)) ==
              ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.CA, p = 0)) ==
              ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.Kt, p = 1)) ==
              ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.GenKt, p = 1.0,
                                     R = 0.4)) == ClusterSequence{PseudoJet}

        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.AntiKt,
                                                   p = 0)
        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.CA,
                                                   p = 1)
        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.Kt,
                                                   p = -1)
        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.GenKt,
                                                   R = 0.4)

        # Supported for now, but will deprecate this calling mode, where only
        # the power is given, at the next major release
        @test typeof(jet_reconstruct(inputs[1]; p = -1)) == ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; p = 0)) == ClusterSequence{PseudoJet}
        @test typeof(jet_reconstruct(inputs[1]; p = 1)) == ClusterSequence{PseudoJet}
        @test_throws KeyError jet_reconstruct(inputs[1]; p = 0.5)

        # No algorithm or power will throw
        @test_throws ArgumentError jet_reconstruct(inputs[1])
    end
end
