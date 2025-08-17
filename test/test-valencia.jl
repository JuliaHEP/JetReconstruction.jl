"""
    test-valencia.jl

Test for the Valencia jet algorithm implementation.

Validates clustering results and cluster sequence against FastJet C++ reference output
using jet matching to handle numerical precision differences.
"""

include("common.jl")

# Reference outputs for R=0.8, beta=1.2, gamma=1.2
valencia_inclusive = ComparisonTestValencia(events_file_ee,
                                            joinpath(@__DIR__, "data",
                                                     "jet-collections-fastjet-valencia-inclusive20-beta1.2-gamma1.2-R0.8-eeH.json.zst"),
                                            JetAlgorithm.Valencia, RecoStrategy.N2Plain,
                                            1.2, 0.8, 1.2,
                                            cs -> inclusive_jets(cs; ptmin = 20.0),
                                            "inclusive beta=1.2 gamma=1.2 R=0.8")
run_reco_test(valencia_inclusive)

valencia_exclusive4 = ComparisonTestValencia(events_file_ee,
                                             joinpath(@__DIR__, "data",
                                                      "jet-collections-fastjet-valencia-exclusive4-beta1.2-gamma1.2-R0.8-eeH.json.zst"),
                                             JetAlgorithm.Valencia, RecoStrategy.N2Plain,
                                             1.2, 0.8, 1.2,
                                             cs -> exclusive_jets(cs; njets = 4),
                                             "exclusive njets beta=1.2 gamma=1.2 R=0.8")
run_reco_test(valencia_exclusive4)

valencia_exclusive_d500 = ComparisonTestValencia(events_file_ee,
                                                 joinpath(@__DIR__, "data",
                                                          "jet-collections-fastjet-valencia-exclusive-d500-beta1.2-gamma1.2-R0.8-eeH.json.zst"),
                                                 JetAlgorithm.Valencia,
                                                 RecoStrategy.N2Plain, 1.2, 0.8, 1.2,
                                                 cs -> exclusive_jets(cs; dcut = 500.0),
                                                 "exclusive dcut=500 beta=1.2 gamma=1.2 R=0.8")
run_reco_test(valencia_exclusive_d500)

#Reference outputs for R=1.0, beta=1.0, gamma=1.0
valencia_inclusive_b1g1 = ComparisonTestValencia(events_file_ee,
                                                 joinpath(@__DIR__, "data",
                                                          "jet-collections-fastjet-valencia-inclusive20-beta1.0-gamma1.0-R1.0-eeH.json.zst"),
                                                 JetAlgorithm.Valencia,
                                                 RecoStrategy.N2Plain, 1.0, 1.0, 1.0,
                                                 cs -> inclusive_jets(cs; ptmin = 20.0),
                                                 "inclusive beta=1.0 gamma=1.0 R=1.0")
run_reco_test(valencia_inclusive_b1g1)

valencia_exclusive4_b1g1 = ComparisonTestValencia(events_file_ee,
                                                  joinpath(@__DIR__, "data",
                                                           "jet-collections-fastjet-valencia-exclusive4-beta1.0-gamma1.0-R1.0-eeH.json.zst"),
                                                  JetAlgorithm.Valencia,
                                                  RecoStrategy.N2Plain, 1.0, 1.0, 1.0,
                                                  cs -> exclusive_jets(cs; njets = 4),
                                                  "exclusive njets beta=1.0 gamma=1.0 R=1.0")
run_reco_test(valencia_exclusive4_b1g1)

valencia_exclusive_d500_b1g1 = ComparisonTestValencia(events_file_ee,
                                                      joinpath(@__DIR__, "data",
                                                               "jet-collections-fastjet-valencia-exclusive-d500-beta1.0-gamma1.0-R1.0-eeH.json.zst"),
                                                      JetAlgorithm.Valencia,
                                                      RecoStrategy.N2Plain, 1.0, 1.0, 1.0,
                                                      cs -> exclusive_jets(cs; dcut = 500.0),
                                                      "exclusive dcut=500 beta=1.0 gamma=1.0 R=1.0")
run_reco_test(valencia_exclusive_d500_b1g1)

# Test dij_dist for Valencia algorithm
@testset "dij_dist Valencia" begin
    # Minimal eereco with two reco jets; only nx,ny,nz,E2p are used by valencia_distance
    eereco = EERecoJet[EERecoJet(1, 0, Inf, Inf, 1.0, 0.0, 0.0, 1.0),
                       EERecoJet(2, 0, Inf, Inf, 0.0, 1.0, 0.0, 2.0)]
    dij = dij_dist(eereco, 1, 2, 1.0, JetAlgorithm.Valencia, 0.8)
    @test dij ≈ valencia_distance(eereco, 1, 2, 0.8)
end

# Valencia distance wrapper coverage
@testset "Valencia distance wrappers" begin
    # Minimal eereco with two reco jets and identical directions so angle=0
    eereco = EERecoJet[EERecoJet(1, 0, Inf, Inf, 1.0, 0.0, 0.0, 9.0),
                       EERecoJet(2, 0, Inf, Inf, 1.0, 0.0, 0.0, 4.0)]
    R = 2.0
    invR2 = inv(R * R)
    @test valencia_distance_inv(eereco, 1, 2, invR2) == 0.0
    @test valencia_distance(eereco, 1, 2, R) == 0.0
end

# Test ee_genkt_algorithm for Valencia algorithm (covers β override and dij_factor selection)
@testset "ee_genkt_algorithm Valencia" begin
    particles = [PseudoJet(1.0, 0.0, 0.0, 1.0)]
    cs = ee_genkt_algorithm(particles; algorithm = JetAlgorithm.Valencia, β = 1.2, γ = 1.2,
                            R = 0.8)
    @test cs isa JetReconstruction.ClusterSequence
end

# Test internal _ee_genkt_algorithm entry (touches StructArray init path too)
@testset "_ee_genkt_algorithm Valencia dij_factor" begin
    # Ensure cluster_hist_index(i) == i for initial jets
    jets = [EEJet(1.0, 0.0, 0.0, 1.0; cluster_hist_index = 1)]
    cs = JetReconstruction._ee_genkt_algorithm(particles = jets,
                                               algorithm = JetAlgorithm.Valencia, p = 1.2,
                                               R = 0.8, γ = 1.2)
    @test cs isa JetReconstruction.ClusterSequence
end


    @testset "Valencia precomputed helpers (unit)" begin
        # Construct a tiny set of EEJet with controlled directions and energies
        # Ensure all jets have non-zero momentum to avoid undefined direction cosines
        p1 = EEJet(0.0, 0.0, 1.0, 10.0; cluster_hist_index = 1)
        p2 = EEJet(0.0, 0.1, 0.5, 5.0; cluster_hist_index = 2)
        p3 = EEJet(0.1, 0.0, -1.0, 2.0; cluster_hist_index = 3)
        parts = [p1, p2, p3]
        N = length(parts)
        R = 1.0
        pwr = 1.2
        γ = 1.0

        # Prepare SoA like Valencia entrypoint
        eereco = StructArray{EERecoJet}(undef, N)
        fill_reco_array!(eereco, parts, R^2, pwr)

        nx = eereco.nx; ny = eereco.ny; nz = eereco.nz; E2p = eereco.E2p
        invR2 = inv(R * R)
        factor = 2 * invR2
        E2p_scaled = similar(E2p)
        beam_term = similar(E2p)
        for k in 1:N
            E2p_scaled[k] = E2p[k] * factor
            sin2 = 1.0 - nz[k] * nz[k]
            beam_term[k] = E2p[k] * sin2^γ
        end

        # valencia_distance_inv_scaled_arrays: compare to manual computation
        i, j = 1, 2
        ang = 1.0 - nx[i]*nx[j] - ny[i]*ny[j] - nz[i]*nz[j]
        @test JetReconstruction.valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j) ≈ min(E2p_scaled[i], E2p_scaled[j]) * ang

        # beam distance arrays
        @test JetReconstruction.valencia_beam_distance_arrays(E2p, nz, 1, γ, pwr) ≈ beam_term[1]

        # initializer: should populate nni/nndist/dijdist without error
        JetReconstruction.get_angular_nearest_neighbours_valencia_precomputed!(eereco, E2p_scaled, beam_term, pwr, γ, R)
        @test all(isfinite, eereco.nndist)
        @test all(isfinite, eereco.dijdist)
        @test all(i -> isa(i, Int), eereco.nni)

        # array-based no-cross and cross update helpers should execute and produce finite outputs
        nnd = copy(eereco.nndist); nni = copy(eereco.nni); dijd = copy(eereco.dijdist)
        JetReconstruction.update_nn_no_cross_arrays_precomputed!(nnd, nni, nx, ny, nz, E2p_scaled, beam_term, dijd, 2, N, 1.0, pwr, γ, R)
        @test isfinite(nnd[2])
        @test isfinite(dijd[2])

        # cross update (should update other slots if appropriate) — assert runs and outputs are finite
        nnd2 = copy(nnd); nni2 = copy(nni); dijd2 = copy(dijd)
        JetReconstruction.update_nn_cross_arrays_precomputed!(nnd2, nni2, nx, ny, nz, E2p_scaled, beam_term, dijd2, 2, N, 1.0, pwr, γ, R)
        @test all(isfinite, nnd2)
        @test all(isfinite, dijd2)
    end
