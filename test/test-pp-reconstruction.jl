# Tests of pp reconstruction algorithms

# Contains all common functions and necessary imports
include("common.jl")

const test_genkt_power = 1.5
const test_cone_size = 0.4

# Test inclusive jets for each algorithm and strategy
for alg in [JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.GenKt],
    stg in [RecoStrategy.N2Plain, RecoStrategy.N2Tiled]

    if alg == JetAlgorithm.GenKt
        power = test_genkt_power
        fastjet_file = joinpath(@__DIR__, "data",
                                "jet-collections-fastjet-inclusive-$(alg)-p$(power).json.zst")
    else
        power = JetReconstruction.algorithm2power[alg]
        fastjet_file = joinpath(@__DIR__, "data",
                                "jet-collections-fastjet-inclusive-$(alg).json.zst")
    end

    test = ComparisonTest(events_file_pp, fastjet_file,
                          alg, stg,
                          power, test_cone_size,
                          (cs) -> inclusive_jets(cs; ptmin = 5.0), "inclusive")
    run_reco_test(test)
end

# Test exclusive njet selections for CA and Kt algorithms
for alg in [JetAlgorithm.CA, JetAlgorithm.Kt]
    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-njets4-$(alg).json.zst")
    test = ComparisonTest(events_file_pp, fastjet_file,
                          alg, RecoStrategy.Best,
                          JetReconstruction.algorithm2power[alg], test_cone_size,
                          (cs) -> exclusive_jets(cs; njets = 4), "exclusive njets")
    run_reco_test(test)
end

# Test exclusive dij selections for CA and Kt algorithms
for (alg, dij_max) in zip([JetAlgorithm.CA, JetAlgorithm.Kt], ["0.99", "20.0"])
    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-dij$(dij_max)-$(alg).json.zst")
    test = ComparisonTest(events_file_pp, fastjet_file,
                          alg, RecoStrategy.Best,
                          JetReconstruction.algorithm2power[alg], test_cone_size,
                          (cs) -> exclusive_jets(cs; dcut = tryparse(Float64, dij_max)),
                          "exclusive dcut")
    run_reco_test(test)
end

# Test for different recombination schemes
# pt-scheme
let
    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-inclusive-AntiKt-ptscheme.json.zst")
    test = ComparisonTest(events_file_pp, fastjet_file,
                          JetAlgorithm.AntiKt, RecoStrategy.N2Tiled,
                          JetReconstruction.algorithm2power[JetAlgorithm.AntiKt],
                          test_cone_size,
                          (cs) -> inclusive_jets(cs; ptmin = 5.0), "AntiKt pt-scheme",
                          addjets_ptscheme, preprocess_ptscheme)
    run_reco_test(test)
end

# pt2-scheme
let
    fastjet_file = joinpath(@__DIR__, "data",
                            "jet-collections-fastjet-inclusive-AntiKt-pt2scheme.json.zst")
    test = ComparisonTest(events_file_pp, fastjet_file,
                          JetAlgorithm.AntiKt, RecoStrategy.N2Tiled,
                          JetReconstruction.algorithm2power[JetAlgorithm.AntiKt],
                          test_cone_size,
                          (cs) -> inclusive_jets(cs; ptmin = 5.0), "AntiKt pt2-scheme",
                          addjets_pt2scheme, preprocess_pt2scheme)
    run_reco_test(test)
end

# Check that all recombination scheme Enums are defined in RecombinationMethods
@testset "RecombinationScheme Enum Check" begin
    for scheme in instances(RecombinationScheme.Recombine)
        @test haskey(RecombinationMethods, scheme)
    end
end
