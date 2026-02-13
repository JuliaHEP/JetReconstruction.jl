# Test the different interface modes of the main jet_reconstruct() function

include("common.jl")

let inputs = JetReconstruction.read_final_state_particles(events_file_ee)
    @testset "jet_reconstruct() interface" begin
        # No algorithm and no power will throw
        @test_throws UndefKeywordError jet_reconstruct(inputs[1])
        # Power but no algorithm will throw
        @test_throws UndefKeywordError jet_reconstruct(inputs[1]; p = 0.5)

        # EE Algorithms
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.Durham)) ==
              ClusterSequence{Float64, EEJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.EEKt, p = -1,
                                     R = 1.0)) == ClusterSequence{Float64, EEJet{Float64}}
        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.EEKt)

        # PP Algorithms
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.AntiKt)) ==
              ClusterSequence{Float64, PseudoJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.CA)) ==
              ClusterSequence{Float64, PseudoJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.Kt)) ==
              ClusterSequence{Float64, PseudoJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.AntiKt, p = -1)) ==
              ClusterSequence{Float64, PseudoJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.CA, p = 0)) ==
              ClusterSequence{Float64, PseudoJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.Kt, p = 1)) ==
              ClusterSequence{Float64, PseudoJet{Float64}}
        @test typeof(jet_reconstruct(inputs[1]; algorithm = JetAlgorithm.GenKt, p = 1.0,
                                     R = 0.4)) == ClusterSequence{Float64, PseudoJet{Float64}}
        @test_throws ArgumentError jet_reconstruct(inputs[1];
                                                   algorithm = JetAlgorithm.GenKt, R = 0.4)
    end
end
