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
    @test valencia_distance(eereco, 1, 2, R) == 0.0
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
