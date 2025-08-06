"""
    test-valencia.jl

Test for the Valencia jet algorithm implementation.

Validates clustering results and cluster sequence against FastJet C++ reference output.
"""

# Include common test utilities
include("common.jl")
import JetReconstruction: pt, rapidity, PseudoJet
 
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
        
        # Apply rapidity filter like in C++ reference: SelectorAbsRapMax(4.0)
        # This removes particles with |rapidity| > 4.0
        filtered_event = filter(p -> abs(rapidity(p)) <= 4.0, event)
        @info "Applied rapidity filter |rap| <= 4.0: $(length(event)) -> $(length(filtered_event)) particles"
        
        # Run Valencia algorithm on first event with β=1.2, γ=0.8
        # These are the same parameters used in the FastJet reference
        clusterseq2 = ee_genkt_algorithm(filtered_event, algorithm=JetAlgorithm.Valencia, p=1.2, γ=0.8)
        
        # Basic sanity checks
        @test length(clusterseq2.jets) >= length(filtered_event)
        @test clusterseq2.algorithm == JetAlgorithm.Valencia
        
        # Test exclusive jets
        exclusive_4jets = exclusive_jets(clusterseq2, njets=4)
        @test length(exclusive_4jets) == 4
        
        # Test inclusive jets (should work even if empty)
        inclusive_jets_result = inclusive_jets(clusterseq2, ptmin=10.0)
        @test isa(inclusive_jets_result, Vector)
        
        # === Detailed comparison with FastJet C++ reference values ===
        
        # Test inclusive jets (pT > 20 GeV) against reference
        inclusive_20gev = inclusive_jets(clusterseq2, ptmin=20.0)
        
        # Load reference data from JSON files
        inclusive_ref_file = joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-inclusive20-eeH.json.zst")
        inclusive_ref_data = read_fastjet_outputs(inclusive_ref_file)[1]["jets"]
        
        @info "Testing inclusive jets (pT > 20 GeV) against FastJet reference:"
        @info "Number of inclusive jets found: $(length(inclusive_20gev))"
        @info "FastJet reference: $(length(inclusive_ref_data)) jets expected"
        
        # Sort by pt (descending) for comparison
        inclusive_sorted = sort(inclusive_20gev, by=jet->pt(PseudoJet(px(jet), py(jet), pz(jet), energy(jet))), rev=true)
        ref_sorted = sort(inclusive_ref_data, by=jet->jet["pt"], rev=true)
        
        @info "Our inclusive jets (pT > 20 GeV, sorted by pT):"
        for (i, jet) in enumerate(inclusive_sorted)
            jet_pj = PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
            jet_pt = pt(jet_pj)
            jet_rap = rapidity(jet_pj)
            @info "  Jet $i: pt = $(round(jet_pt, digits=6)), rap = $(round(jet_rap, digits=7))"
        end
        
        @info "FastJet reference inclusive jets:"
        for (i, jet) in enumerate(ref_sorted)
            @info "  Jet $i: pt = $(jet["pt"]), rap = $(jet["rap"])"
        end
        
        # Test against reference with tolerance
        @test length(inclusive_20gev) == length(inclusive_ref_data)
        for (i, (our_jet, ref_jet)) in enumerate(zip(inclusive_sorted, ref_sorted))
            jet_pj = PseudoJet(px(our_jet), py(our_jet), pz(our_jet), energy(our_jet))
            our_pt = pt(jet_pj)
            our_rap = rapidity(jet_pj)
            
            @test abs(our_pt - ref_jet["pt"]) < 0.1  # 0.1 GeV tolerance
            @test abs(our_rap - ref_jet["rap"]) < 0.01  # 0.01 rapidity tolerance
        end
        
        # Test exclusive N=4 clustering against reference
        exclusive_4 = exclusive_jets(clusterseq2, njets=4)
        @test length(exclusive_4) == 4
        
        # Load reference data from JSON file
        exclusive_ref_file = joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-exclusive4-eeH.json.zst")
        exclusive_ref_data = read_fastjet_outputs(exclusive_ref_file)[1]["jets"]
        
        # Convert to PseudoJets for pt/rapidity calculation and sort by pt (descending)
        jets_with_pt = [(PseudoJet(px(jet), py(jet), pz(jet), energy(jet)), i) for (i, jet) in enumerate(exclusive_4)]
        sort!(jets_with_pt, by=x->pt(x[1]), rev=true)
        ref_sorted = sort(exclusive_ref_data, by=jet->jet["pt"], rev=true)
        
        # Print our results for comparison
        @info "Valencia algorithm comparison with FastJet reference (exclusive N=4, with |rap| <= 4.0 filter):"
        @info "Our Valencia results (exclusive N=4, sorted by pt):"
        for i in 1:4
            jet_pt = pt(jets_with_pt[i][1])
            jet_rap = rapidity(jets_with_pt[i][1])
            @info "  Jet $i: pt = $(round(jet_pt, digits=6)), rap = $(round(jet_rap, digits=7))"
        end
        
        @info "FastJet reference (exclusive N=4):"
        for (i, jet) in enumerate(ref_sorted)
            @info "  Jet $i: pt = $(jet["pt"]), rap = $(jet["rap"])"
        end
        
        # Test against reference with tolerance
        @test length(exclusive_4) == length(exclusive_ref_data)
        for (i, (our_jet_data, ref_jet)) in enumerate(zip(jets_with_pt, ref_sorted))
            our_pt = pt(our_jet_data[1])
            our_rap = rapidity(our_jet_data[1])
            
            @test abs(our_pt - ref_jet["pt"]) < 0.1  # 0.1 GeV tolerance
            @test abs(our_rap - ref_jet["rap"]) < 0.01  # 0.01 rapidity tolerance
        end
        
        # Test exclusive clustering up to d = 500 against reference
        exclusive_d500 = exclusive_jets(clusterseq2, dcut=500.0)
        
        # Load reference data from JSON file
        exclusive_d500_ref_file = joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-exclusive-d500-eeH.json.zst")
        exclusive_d500_ref_data = read_fastjet_outputs(exclusive_d500_ref_file)[1]["jets"]
        
        @info "Testing exclusive clustering up to d = 500 against FastJet reference:"
        @info "Number of exclusive jets found: $(length(exclusive_d500))"
        @info "FastJet reference: $(length(exclusive_d500_ref_data)) jets expected"
        
        # Sort by pt (descending) for comparison
        exclusive_d500_sorted = sort(exclusive_d500, by=jet->pt(PseudoJet(px(jet), py(jet), pz(jet), energy(jet))), rev=true)
        ref_d500_sorted = sort(exclusive_d500_ref_data, by=jet->jet["pt"], rev=true)
        
        @info "Our exclusive jets (d < 500, sorted by pT):"
        for (i, jet) in enumerate(exclusive_d500_sorted)
            jet_pj = PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
            jet_pt = pt(jet_pj)
            jet_rap = rapidity(jet_pj)
            @info "  Jet $i: pt = $(round(jet_pt, digits=6)), rap = $(round(jet_rap, digits=7))"
        end
        
        @info "FastJet reference exclusive jets (d < 500):"
        for (i, jet) in enumerate(ref_d500_sorted)
            @info "  Jet $i: pt = $(jet["pt"]), rap = $(jet["rap"])"
        end
        
        # Test against reference with tolerance
        @test length(exclusive_d500) == length(exclusive_d500_ref_data)
        for (i, (our_jet, ref_jet)) in enumerate(zip(exclusive_d500_sorted, ref_d500_sorted))
            jet_pj = PseudoJet(px(our_jet), py(our_jet), pz(our_jet), energy(our_jet))
            our_pt = pt(jet_pj)
            our_rap = rapidity(jet_pj)

            @test abs(our_pt - ref_jet["pt"]) < 0.1  # 0.1 GeV tolerance
            @test abs(our_rap - ref_jet["rap"]) < 0.01  # 0.01 rapidity tolerance
        end
        
    end
end
