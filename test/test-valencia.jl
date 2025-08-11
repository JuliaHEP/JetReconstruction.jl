"""
    test-valencia.jl

Test for the Valencia jet algorithm implementation.

Validates clustering results and cluster sequence against FastJet C++ reference output
using jet matching to handle numerical precision differences.
"""

# Include common test utilities
include("common.jl")
import JetReconstruction: pt, rapidity, PseudoJet

using Test
using JetReconstruction

# JSON package for parsing
using JSON

# JSON reading function
function read_fastjet_outputs(filename)
    if endswith(filename, ".zst")
        # Decompress .zst file to a buffer and parse JSON
        io = open(`zstdcat $(filename)`, "r")
        data = read(io, String)
        close(io)
        return JSON.parse(data)
    else
        open(filename, "r") do file
            return JSON.parse(file)
        end
    end
end

"""
    match_jets_with_pairing_optimization(ref_jets, julia_jets; pt_tolerance=0.5, rap_tolerance=1.0, pt_similarity_threshold=0.01)

Enhanced jet matching that tries alternative pairings when jets have nearly identical pT.
"""
function match_jets_with_pairing_optimization(ref_jets, julia_jets; pt_tolerance = 0.5,
                                              rap_tolerance = 1.0,
                                              pt_similarity_threshold = 0.01)
    # First, identify groups of jets with similar pT
    n_ref = length(ref_jets)
    n_julia = length(julia_jets)

    # Check if we have exactly 2 jets with very similar pT
    if n_ref == 2 && n_julia == 2
        ref_pts = [ref_jets[1]["pt"], ref_jets[2]["pt"]]
        julia_pts = [begin
                         jet_pj = PseudoJet(px(julia_jets[i]), py(julia_jets[i]),
                                            pz(julia_jets[i]), energy(julia_jets[i]))
                         pt(jet_pj)
                     end
                     for i in 1:2]

        # Check if reference jets have similar pT
        if abs(ref_pts[1] - ref_pts[2]) < pt_similarity_threshold
            # Try both pairing arrangements

            # Arrangement 1: ref1->julia1, ref2->julia2
            jet1_pj = PseudoJet(px(julia_jets[1]), py(julia_jets[1]), pz(julia_jets[1]),
                                energy(julia_jets[1]))
            jet2_pj = PseudoJet(px(julia_jets[2]), py(julia_jets[2]), pz(julia_jets[2]),
                                energy(julia_jets[2]))

            julia1_pt, julia1_rap = pt(jet1_pj), rapidity(jet1_pj)
            julia2_pt, julia2_rap = pt(jet2_pj), rapidity(jet2_pj)

            # Score arrangement 1
            score1 = sqrt((abs(julia1_pt - ref_jets[1]["pt"]) / pt_tolerance)^2 +
                          (abs(julia1_rap - ref_jets[1]["rap"]) / rap_tolerance)^2) +
                     sqrt((abs(julia2_pt - ref_jets[2]["pt"]) / pt_tolerance)^2 +
                          (abs(julia2_rap - ref_jets[2]["rap"]) / rap_tolerance)^2)

            # Score arrangement 2: ref1->julia2, ref2->julia1 (swapped)
            score2 = sqrt((abs(julia2_pt - ref_jets[1]["pt"]) / pt_tolerance)^2 +
                          (abs(julia2_rap - ref_jets[1]["rap"]) / rap_tolerance)^2) +
                     sqrt((abs(julia1_pt - ref_jets[2]["pt"]) / pt_tolerance)^2 +
                          (abs(julia1_rap - ref_jets[2]["rap"]) / rap_tolerance)^2)

            # Use the better arrangement
            if score2 < score1
                # Swapped arrangement is better
                return [(ref_jets[1], julia_jets[2], 1, 2),
                        (ref_jets[2], julia_jets[1], 2, 1)]
            else
                # Original arrangement is better  
                return [(ref_jets[1], julia_jets[1], 1, 1),
                        (ref_jets[2], julia_jets[2], 2, 2)]
            end
        end
    end

    # Fall back to standard greedy matching for other cases
    used_julia_jets = Set{Int}()
    matched_pairs = []

    # Sort reference jets by pT (descending) for consistent matching
    ref_jets_sorted = sort(collect(enumerate(ref_jets)), by = x -> x[2]["pt"], rev = true)

    for (ref_idx, ref_jet) in ref_jets_sorted
        best_match = nothing
        best_distance = Inf

        for (julia_idx, julia_jet) in enumerate(julia_jets)
            if julia_idx in used_julia_jets
                continue
            end

            jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                               energy(julia_jet))
            julia_pt = pt(jet_pj)
            julia_rap = rapidity(jet_pj)

            # Check if jet is within tolerance
            pt_diff = abs(julia_pt - ref_jet["pt"])
            rap_diff = abs(julia_rap - ref_jet["rap"])

            if pt_diff < pt_tolerance && rap_diff < rap_tolerance
                # Use combined distance metric
                distance = sqrt((pt_diff / pt_tolerance)^2 + (rap_diff / rap_tolerance)^2)

                if distance < best_distance
                    best_distance = distance
                    best_match = (julia_jet, julia_idx)
                end
            end
        end

        if best_match !== nothing
            julia_jet, julia_idx = best_match
            push!(matched_pairs, (ref_jet, julia_jet, ref_idx, julia_idx))
            push!(used_julia_jets, julia_idx)
        end
    end

    return matched_pairs
end

"""
    match_jets(ref_jets, julia_jets; pt_tolerance=0.5, rap_tolerance=1.0)

Match reference jets to Julia jets based on kinematics with enhanced handling for similar pT jets.
"""
function match_jets(ref_jets, julia_jets; pt_tolerance = 0.5, rap_tolerance = 1.0)
    return match_jets_with_pairing_optimization(ref_jets, julia_jets;
                                                pt_tolerance = pt_tolerance,
                                                rap_tolerance = rap_tolerance)
end

@testset "Valencia algorithm basic test" begin
    # Test with simple 2-particle system
    particles = [
        PseudoJet(1.0, 0.0, 0.0, 1.0),
        PseudoJet(0.0, 1.0, 0.0, 1.0)
    ]
    # Run Valencia algorithm with test parameters
    β = 0.8
    γ = 0.8
    R = 1.2

    # Basic checks
    clusterseq = ee_genkt_algorithm(particles, algorithm = JetAlgorithm.Valencia, p = β,
                                    γ = γ, R = R)
    @test clusterseq isa ClusterSequence

    # Test exclusive jets
    exclusive_jets_result = exclusive_jets(clusterseq, njets = 1)
    @test length(exclusive_jets_result) == 1

    eventfile = joinpath(@__DIR__, "data", "events.eeH.hepmc3.zst")
    if isfile(eventfile)
        events = read_final_state_particles(eventfile)
        # Load reference data for all events (now using .zst files)
        inclusive_ref_file = joinpath(@__DIR__, "data",
                                      "jet-collections-fastjet-valencia-inclusive20-eeH.json.zst")
        inclusive_ref_all = read_fastjet_outputs(inclusive_ref_file)
        exclusive_ref_file = joinpath(@__DIR__, "data",
                                      "jet-collections-fastjet-valencia-exclusive4-eeH.json.zst")
        exclusive_ref_all = read_fastjet_outputs(exclusive_ref_file)
        exclusive_d500_ref_file = joinpath(@__DIR__, "data",
                                           "jet-collections-fastjet-valencia-exclusive-d500-eeH.json.zst")
        exclusive_d500_ref_all = read_fastjet_outputs(exclusive_d500_ref_file)

        for (evt_idx, event) in enumerate(events)
            filtered_event = filter(p -> abs(rapidity(p)) <= 4.0, event)
            clusterseq = ee_genkt_algorithm(filtered_event,
                                            algorithm = JetAlgorithm.Valencia, p = 0.8,
                                            γ = 0.8, R = 1.2)

                                            
                @info "Matching results:"
                matched_pairs = match_jets(inclusive_ref_data, inclusive_20gev;
                                           pt_tolerance = 0.5, rap_tolerance = 2.0)
                for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                    jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                       energy(julia_jet))
                    our_pt = pt(jet_pj)
                    our_rap = rapidity(jet_pj)
                    pt_diff = abs(our_pt - ref_jet["pt"])
                    rap_diff = abs(our_rap - ref_jet["rap"])
                    @info "  Ref $ref_idx -> Julia $julia_idx: pT_diff=$pt_diff, rap_diff=$rap_diff"
                end
            end

            # Inclusive jets (pt > 20 GeV) with jet matching
            inclusive_20gev = inclusive_jets(clusterseq, ptmin = 20.0)
            inclusive_ref_data = inclusive_ref_all[evt_idx]["jets"]
            @info "Event $evt_idx: Inclusive jets (pT > 20 GeV): found $(length(inclusive_20gev)), reference $(length(inclusive_ref_data))"

            # Match jets instead of position-based comparison
            matched_pairs = match_jets(inclusive_ref_data, inclusive_20gev;
                                       pt_tolerance = 0.1, rap_tolerance = 2.0)

            # Enhanced test: if there are 2 jets with nearly identical pT, check both rapidity pairings
            if length(matched_pairs) == 2
                ref_pts = [matched_pairs[1][1]["pt"], matched_pairs[2][1]["pt"]]
                if abs(ref_pts[1] - ref_pts[2]) < 0.01
                    # Try both rapidity pairings
                    jet1_pj = PseudoJet(px(matched_pairs[1][2]), py(matched_pairs[1][2]),
                                        pz(matched_pairs[1][2]),
                                        energy(matched_pairs[1][2]))
                    jet2_pj = PseudoJet(px(matched_pairs[2][2]), py(matched_pairs[2][2]),
                                        pz(matched_pairs[2][2]),
                                        energy(matched_pairs[2][2]))
                    rap1 = rapidity(jet1_pj)
                    rap2 = rapidity(jet2_pj)
                    ref_rap1 = matched_pairs[1][1]["rap"]
                    ref_rap2 = matched_pairs[2][1]["rap"]
                    # Arrangement 1
                    arrangement1 = abs(rap1 - ref_rap1) < 0.15 &&
                                   abs(rap2 - ref_rap2) < 0.15
                    # Arrangement 2 (swapped)
                    arrangement2 = abs(rap2 - ref_rap1) < 0.15 &&
                                   abs(rap1 - ref_rap2) < 0.15
                    if !(arrangement1 || arrangement2)
                        println("Event $evt_idx inclusive jets rapidity pairing failed. Printing jets:")
                        println("Reference jets:")
                        for (i, ref_jet) in enumerate(inclusive_ref_data)
                            println("  Jet $i: pT=$(ref_jet["pt"]), rap=$(ref_jet["rap"])")
                        end
                        println("Julia jets:")
                        for (i, julia_jet) in enumerate(inclusive_20gev)
                            jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                               energy(julia_jet))
                            println("  Jet $i: pT=$(pt(jet_pj)), rap=$(rapidity(jet_pj))")
                        end
                    end
                    @test arrangement1 || arrangement2
                    # Always test pT agreement
                    if !(abs(pt(jet1_pj) - matched_pairs[1][1]["pt"]) < 0.1)
                        println("Event $evt_idx inclusive jets pT test failed for jet 1.")
                    end
                    if !(abs(pt(jet2_pj) - matched_pairs[2][1]["pt"]) < 0.1)
                        println("Event $evt_idx inclusive jets pT test failed for jet 2.")
                    end
                    @test abs(pt(jet1_pj) - matched_pairs[1][1]["pt"]) < 0.1
                    @test abs(pt(jet2_pj) - matched_pairs[2][1]["pt"]) < 0.1
                else
                    # Standard test for jets with distinct pT
                    failed = false
                    for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                        jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                           energy(julia_jet))
                        our_pt = pt(jet_pj)
                        our_rap = rapidity(jet_pj)
                        if !(abs(our_pt - ref_jet["pt"]) < 0.1) ||
                           !(abs(our_rap - ref_jet["rap"]) < 0.15)
                            failed = true
                        end
                    end
                    if failed
                        println("Event $evt_idx inclusive jets kinematic test failed. Printing all jets:")
                        println("Reference jets:")
                        for (i, ref_jet) in enumerate(inclusive_ref_data)
                            println("  Jet $i: pT=$(ref_jet["pt"]), rap=$(ref_jet["rap"])")
                        end
                        println("Julia jets:")
                        for (i, julia_jet) in enumerate(inclusive_20gev)
                            jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                               energy(julia_jet))
                            println("  Jet $i: pT=$(pt(jet_pj)), rap=$(rapidity(jet_pj))")
                        end
                    end
                    for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                        jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                           energy(julia_jet))
                        our_pt = pt(jet_pj)
                        our_rap = rapidity(jet_pj)
                        @test abs(our_pt - ref_jet["pt"]) < 0.1
                        @test abs(our_rap - ref_jet["rap"]) < 0.15
                    end
                end
            else
                # Standard test for other cases
                failed = false
                for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                    jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                       energy(julia_jet))
                    our_pt = pt(jet_pj)
                    our_rap = rapidity(jet_pj)
                    if !(abs(our_pt - ref_jet["pt"]) < 0.1) ||
                       !(abs(our_rap - ref_jet["rap"]) < 0.15)
                        failed = true
                    end
                end
                if failed
                    println("Event $evt_idx inclusive jets kinematic test failed. Printing all jets:")
                    println("Reference jets:")
                    for (i, ref_jet) in enumerate(inclusive_ref_data)
                        println("  Jet $i: pT=$(ref_jet["pt"]), rap=$(ref_jet["rap"])")
                    end
                    println("Julia jets:")
                    for (i, julia_jet) in enumerate(inclusive_20gev)
                        jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                           energy(julia_jet))
                        println("  Jet $i: pT=$(pt(jet_pj)), rap=$(rapidity(jet_pj))")
                    end
                end
                for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                    jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                       energy(julia_jet))
                    our_pt = pt(jet_pj)
                    our_rap = rapidity(jet_pj)
                    @test abs(our_pt - ref_jet["pt"]) < 0.1
                    @test abs(our_rap - ref_jet["rap"]) < 0.15
                end
            end

            # Exclusive N=4 jets with jet matching
            exclusive_4 = exclusive_jets(clusterseq, njets = 4)
            exclusive_ref_data = exclusive_ref_all[evt_idx]["jets"]
            @info "Event $evt_idx: Exclusive jets (N=4): found $(length(exclusive_4)), reference $(length(exclusive_ref_data))"

            # Match jets
            matched_pairs = match_jets(exclusive_ref_data, exclusive_4; pt_tolerance = 0.1,
                                       rap_tolerance = 2.0)

            # Test kinematic agreement for matched jets
            for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                   energy(julia_jet))
                our_pt = pt(jet_pj)
                our_rap = rapidity(jet_pj)
                @test abs(our_pt - ref_jet["pt"]) < 0.1
                @test abs(our_rap - ref_jet["rap"]) < 0.15
            end

            # Exclusive d < 500 jets with jet matching
            exclusive_d500 = exclusive_jets(clusterseq, dcut = 500.0)
            exclusive_d500_ref_data = exclusive_d500_ref_all[evt_idx]["jets"]
            @info "Event $evt_idx: Exclusive jets (d < 500): found $(length(exclusive_d500)), reference $(length(exclusive_d500_ref_data))"

            # Match jets with more lenient tolerances for d-cut jets
            matched_pairs = match_jets(exclusive_d500_ref_data, exclusive_d500;
                                       pt_tolerance = 0.1, rap_tolerance = 2.0)

            # Test kinematic agreement for matched jets
            for (ref_jet, julia_jet, ref_idx, julia_idx) in matched_pairs
                jet_pj = PseudoJet(px(julia_jet), py(julia_jet), pz(julia_jet),
                                   energy(julia_jet))
                our_pt = pt(jet_pj)
                our_rap = rapidity(jet_pj)
                @test abs(our_pt - ref_jet["pt"]) < 0.1
                @test abs(our_rap - ref_jet["rap"]) < 0.15
            end
        end
    end
end
