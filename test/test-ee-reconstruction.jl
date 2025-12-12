# Tests of e+e- reconstruction algorithms

include("common.jl")

durham_njets3 = ComparisonTest(events_file_ee,
                               joinpath(@__DIR__, "data",
                                        "jet-collections-fastjet-njets3-Durham-eeH.json.zst"),
                               JetAlgorithm.Durham, RecoStrategy.N2Plain, 1, 4.0,
                               (cs) -> exclusive_jets(cs; njets = 3), "exclusive njets")

run_reco_test(durham_njets3)

eekt_inclusive = ComparisonTest(events_file_ee,
                                joinpath(@__DIR__, "data",
                                         "jet-collections-fastjet-inclusive-EEKt-p-1-R1.0.json.zst"),
                                JetAlgorithm.EEKt, RecoStrategy.N2Plain, -1, 1.0,
                                (cs) -> inclusive_jets(cs; ptmin = 5.0), "inclusive")

run_reco_test(eekt_inclusive)

for r in [2.0, 4.0]
    eekt_njets = ComparisonTest(events_file_ee,
                                joinpath(@__DIR__, "data",
                                         "jet-collections-fastjet-njets4-EEKt-p1-R$(r).json.zst"),
                                JetAlgorithm.EEKt, RecoStrategy.N2Plain, 1, r,
                                (cs) -> exclusive_jets(cs; njets = 4), "exclusive njets")
    run_reco_test(eekt_njets)
end

# Optimization/helper coverage

@testset "E2p integer fast paths in fill/insert (e+e-)" begin
    E1, E2 = 3.0, 2.0
    # PseudoJet signature is (px, py, pz, E); give nonzero momentum to avoid 1/0 in EEJet
    pj1 = PseudoJet(1.0, 0.0, 0.0, E1)
    pj2 = PseudoJet(1.0, 0.0, 0.0, E2)
    particles = EEJet[]
    push!(particles, EEJet(pj1; cluster_hist_index = 1))
    push!(particles, EEJet(pj2; cluster_hist_index = 2))

    eereco = StructArray{EERecoJet}(undef, 2)
    R = 1.0
    R2 = R^2

    fill_reco_array!(eereco, particles, R2, 1)
    @test eereco.E2p[1] ≈ E1^2
    @test eereco.E2p[2] ≈ E2^2

    fill_reco_array!(eereco, particles, R2, 2)
    @test eereco.E2p[1] ≈ E1^4
    @test eereco.E2p[2] ≈ E2^4

    merged = EEJet(1.0, 0.0, 0.0, 5.0; cluster_hist_index = 0)
    insert_new_jet!(eereco, 1, 3, R2, merged, 1)
    @test eereco.E2p[1] ≈ 5.0^2

    insert_new_jet!(eereco, 2, 4, R2, merged, 2)
    @test eereco.E2p[2] ≈ 5.0^4
end

@testset "copy_to_slot! copies fields" begin
    eereco = StructArray{EERecoJet}(undef, 2)
    eereco.index[1] = 10
    eereco.nni[1] = 1
    eereco.nndist[1] = 3.14
    eereco.dijdist[1] = 2.71
    eereco.nx[1] = 0.1
    eereco.ny[1] = 0.2
    eereco.nz[1] = 0.3
    eereco.E2p[1] = 7.0
    copy_to_slot!(eereco, 1, 2)
    @test eereco.index[2] == 10
    @test eereco.dijdist[2] == 2.71
    @test eereco.E2p[2] == 7.0
end
