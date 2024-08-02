# Tests of e+e- reconstruction algorithms

include("common.jl")

durham_njets3 = ComparisonTest(events_file_ee,
                               joinpath(@__DIR__, "data", "jet-collections-fastjet-njets3-Durham-eeH.json.gz"),
                               JetAlgorithm.Durham, RecoStrategy.N2Plain, 1, 4.0,
                               (cs) -> exclusive_jets(cs; njets = 3))
                               
run_reco_test(durham_njets3)
