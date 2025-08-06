"""
Test for Valencia algorithm implementation
"""

using Test
using JetReconstruction

@testset "Valencia algorithm basic test" begin
    # Test with simple 2-particle system
    particles = [
        PseudoJet(1.0, 0.0, 0.0, 1.0),
        PseudoJet(0.0, 1.0, 0.0, 1.0)
    ]
    
    # Run Valencia algorithm with test parameters
    β = 1.2
    γ = 0.8
    clusterseq = ee_genkt_algorithm(particles, algorithm=JetAlgorithm.Valencia, p=β, γ=γ)
    
    # Basic checks
    @test length(clusterseq.jets) >= 2  # At least input particles
    @test clusterseq.algorithm == JetAlgorithm.Valencia
    @test clusterseq.power == β
    
    # Test exclusive jets
    exclusive_jets_result = exclusive_jets(clusterseq, njets=1)
    @test length(exclusive_jets_result) == 1
    
    # Test the algorithm can run without errors on real data
    eventfile = joinpath(@__DIR__, "data", "events.eeH.hepmc3.zst")
    if isfile(eventfile)
        events = read_final_state_particles(eventfile)
        event = events[1]
        
        # Run Valencia algorithm on first event
        clusterseq2 = ee_genkt_algorithm(event, algorithm=JetAlgorithm.Valencia, p=1.2, γ=0.8)
        
        # Basic sanity checks
        @test length(clusterseq2.jets) >= length(event)
        @test clusterseq2.algorithm == JetAlgorithm.Valencia
        
        # Test exclusive jets
        exclusive_4jets = exclusive_jets(clusterseq2, njets=4)
        @test length(exclusive_4jets) == 4
        
        # Test inclusive jets (should work even if empty)
        inclusive_jets_result = inclusive_jets(clusterseq2, ptmin=10.0)
        @test isa(inclusive_jets_result, Vector)
    end
end
