# Tests of the winner-takes-all (WTA) recombination schemes

include("common.jl")

# Comparisons with FastJet's WTA_pt_scheme and WTA_modp_scheme on the same events
for (scheme, label) in ((RecombinationScheme.WTAPtScheme, "wtaptscheme"),
                        (RecombinationScheme.WTAModpScheme, "wtamodpscheme"))
    methods = RecombinationMethods[scheme]
    for strategy in [RecoStrategy.N2Plain, RecoStrategy.N2Tiled]
        fastjet_file = joinpath(@__DIR__, "data",
                                "jet-collections-fastjet-inclusive-AntiKt-$(label).json.zst")
        test = ComparisonTest(events_file_pp, fastjet_file,
                              JetAlgorithm.AntiKt, strategy,
                              JetReconstruction.algorithm2power[JetAlgorithm.AntiKt], 0.4,
                              (cs) -> inclusive_jets(cs; ptmin = 5.0), "AntiKt $(scheme)",
                              methods.recombine, methods.preprocess)
        run_reco_test(test)
    end

    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-njets3-Durham-$(label)-eeH.json.zst")
    test = ComparisonTest(events_file_ee, fastjet_file,
                          JetAlgorithm.Durham, RecoStrategy.N2Plain, 1, 4.0,
                          (cs) -> exclusive_jets(cs; njets = 3),
                          "exclusive njets $(scheme)",
                          methods.recombine, methods.preprocess)
    run_reco_test(test)
end

# Unit vector along the 3-momentum of a jet
function direction(jet)
    [JetReconstruction.px(jet), JetReconstruction.py(jet),
        JetReconstruction.pz(jet)] ./ JetReconstruction.p(jet)
end

@testset "WTA recombination properties ($(T))" for T in (PseudoJet, EEJet)
    # `hard` has both the larger pt and the larger |p|; both are massive
    hard = T(1.0, 2.0, 30.0, 40.0)
    soft = T(-0.5, 0.2, 3.0, 3.1)

    for recombine in (addjets_wta_ptscheme, addjets_wta_modpscheme)
        merged = @inferred recombine(hard, soft)
        @test merged isa T
        @test JetReconstruction.cluster_hist_index(recombine(hard, soft;
                                                             cluster_hist_index = 7)) == 7
        # The harder jet wins whatever the argument order, and gives its mass
        swapped = recombine(soft, hard)
        @test JetReconstruction.energy(swapped) ≈ JetReconstruction.energy(merged)
        @test direction(swapped) ≈ direction(merged)
        @test JetReconstruction.mass(merged) ≈ JetReconstruction.mass(hard)
    end

    merged = addjets_wta_ptscheme(hard, soft)
    @test JetReconstruction.rapidity(merged) ≈ JetReconstruction.rapidity(hard)
    @test JetReconstruction.phi(merged) ≈ JetReconstruction.phi(hard)
    @test JetReconstruction.pt(merged) ≈
          JetReconstruction.pt(hard) + JetReconstruction.pt(soft)

    merged = addjets_wta_modpscheme(hard, soft)
    @test direction(merged) ≈ direction(hard)
    @test JetReconstruction.p(merged) ≈
          JetReconstruction.p(hard) + JetReconstruction.p(soft)

    # On a tie in the ordering variable the first jet wins
    massive = T(3.0, 0.0, 0.0, 5.0)    # pt = |p| = 3, m = 4
    massless = T(0.0, 3.0, 0.0, 3.0)   # pt = |p| = 3, m = 0
    for recombine in (addjets_wta_ptscheme, addjets_wta_modpscheme)
        @test JetReconstruction.mass(recombine(massive, massless)) ≈ 4.0
        @test JetReconstruction.mass2(recombine(massless, massive))≈0.0 atol=1e-9
    end

    # If the harder jet is at rest, so is the merged jet, with E = m
    at_rest = addjets_wta_modpscheme(T(0.0, 0.0, 0.0, 2.0), T(0.0, 0.0, 0.0, 1.0))
    @test (JetReconstruction.px(at_rest), JetReconstruction.py(at_rest),
           JetReconstruction.pz(at_rest), JetReconstruction.energy(at_rest)) ==
          (0.0, 0.0, 0.0, 2.0)
end
