# Tests of reconstruction algorithms in Float32

include("common.jl")

durham_njets3_f32 = ComparisonTest(events_file_ee,
                                   joinpath(@__DIR__, "data",
                                            "jet-collections-fastjet-njets3-Durham-eeH.json.gz"),
                                   JetAlgorithm.Durham, RecoStrategy.N2Plain, 1, 4.0,
                                   (cs) -> exclusive_jets(cs; njets = 3),
                                   "exclusive njets Float32", Float32)

run_reco_test(durham_njets3_f32)

antikt_pcut_f32 = ComparisonTest(events_file_pp,
                                 joinpath(@__DIR__, "data",
                                          "jet-collections-fastjet-inclusive-AntiKt.json.gz"),
                                 JetAlgorithm.AntiKt, RecoStrategy.N2Tiled, -1, 0.4,
                                 (cs) -> inclusive_jets(cs; ptmin = 5.0),
                                 "inclusive-float32", Float32)
run_reco_test(antikt_pcut_f32)
