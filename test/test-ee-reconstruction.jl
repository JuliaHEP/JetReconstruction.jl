# Tests of e+e- reconstruction algorithms

include("common.jl")

durham_njets3 = ComparisonTest(events_file_ee,
                               joinpath(@__DIR__, "data",
                                        "jet-collections-fastjet-njets3-Durham-eeH.json.gz"),
                               JetAlgorithm.Durham, RecoStrategy.N2Plain, 1, 4.0,
                               (cs) -> exclusive_jets(cs; njets = 3), "exclusive njets")

run_reco_test(durham_njets3)

eekt_inclusive = ComparisonTest(events_file_ee,
                                joinpath(@__DIR__, "data",
                                         "jet-collections-fastjet-inclusive-EEKt-p-1-R1.0.json.gz"),
                                JetAlgorithm.EEKt, RecoStrategy.N2Plain, -1, 1.0,
                                (cs) -> inclusive_jets(cs; ptmin = 5.0), "inclusive")

run_reco_test(eekt_inclusive)

for r in [2.0, 4.0]
    eekt_njets = ComparisonTest(events_file_ee,
                                joinpath(@__DIR__, "data",
                                         "jet-collections-fastjet-njets4-EEKt-p1-R$(r).json.gz"),
                                JetAlgorithm.EEKt, RecoStrategy.N2Plain, 1, r,
                                (cs) -> exclusive_jets(cs; njets = 4), "exclusive njets")
    run_reco_test(eekt_njets)
end
