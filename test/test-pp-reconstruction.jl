# Tests of pp reconstruction algorithms

# Contains all common functions and necessary imports
include("common.jl")

# Test inclusive jets for each algorithm and strategy
for alg in [JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt],
    stg in [RecoStrategy.N2Plain, RecoStrategy.N2Tiled]

    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-inclusive-$(alg).json.gz")
    test = ComparisonTest(events_file_pp, fastjet_file,
                        alg, stg,
                        JetReconstruction.algorithm2power[alg], 0.4,
                        (cs) -> inclusive_jets(cs; ptmin = 5.0))
    run_reco_test(test)
end
