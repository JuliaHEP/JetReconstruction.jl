# Tests of pp reconstruction algorithms

# Contains all common functions and necessary imports
include("common.jl")

# Test inclusive jets for each algorithm and strategy
for alg in [JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.GenKt],
    stg in [RecoStrategy.N2Plain, RecoStrategy.N2Tiled]

    if alg == JetAlgorithm.GenKt
        power = 1.5
        fastjet_file = joinpath(@__DIR__, "data",
                                "jet-collections-fastjet-inclusive-$(alg)-p$(power).json.gz")
    else
        power = JetReconstruction.algorithm2power[alg]
        fastjet_file = joinpath(@__DIR__, "data",
                                "jet-collections-fastjet-inclusive-$(alg).json.gz")
    end

    test = ComparisonTest(events_file_pp, fastjet_file,
                          alg, stg,
                          power, 0.4,
                          (cs) -> inclusive_jets(cs; ptmin = 5.0), "inclusive")
    run_reco_test(test)
end

# Test exclusive njet selections for CA and Kt algorithms
for alg in [JetAlgorithm.CA, JetAlgorithm.Kt]
    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-njets4-$(alg).json.gz")
    test = ComparisonTest(events_file_pp, fastjet_file,
                          alg, RecoStrategy.Best,
                          JetReconstruction.algorithm2power[alg], 0.4,
                          (cs) -> exclusive_jets(cs; njets = 4), "exclusive njets")
    run_reco_test(test)
end

# Test exclusive dij selections for CA and Kt algorithms
for (alg, dij_max) in zip([JetAlgorithm.CA, JetAlgorithm.Kt], ["0.99", "20.0"])
    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-dij$(dij_max)-$(alg).json.gz")
    test = ComparisonTest(events_file_pp, fastjet_file,
                          alg, RecoStrategy.Best,
                          JetReconstruction.algorithm2power[alg], 0.4,
                          (cs) -> exclusive_jets(cs; dcut = tryparse(Float64, dij_max)),
                          "exclusive dcut")
    run_reco_test(test)
end
