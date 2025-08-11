"""
    test-valencia.jl

Test for the Valencia jet algorithm implementation.

Validates clustering results and cluster sequence against FastJet C++ reference output
using jet matching to handle numerical precision differences.
"""

# Include common test utilities
include("common.jl")
import JetReconstruction: pt, rapidity, PseudoJet

using Test
using JetReconstruction

# Reference outputs for R=0.8, beta=1.2, gamma=1.2
valencia_inclusive = ComparisonTestValencia(
    events_file_ee,
    joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-inclusive20-beta1.2-gamma1.2-R0.8-eeH.json.zst"),
    JetAlgorithm.Valencia, RecoStrategy.N2Plain, 1.2, 0.8, 1.2,
    cs -> inclusive_jets(cs; ptmin = 20.0), "inclusive beta=1.2 gamma=1.2 R=0.8"
)
run_reco_test(valencia_inclusive)

valencia_exclusive4 = ComparisonTestValencia(
    events_file_ee,
    joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-exclusive4-beta1.2-gamma1.2-R0.8-eeH.json.zst"),
    JetAlgorithm.Valencia, RecoStrategy.N2Plain, 1.2, 0.8, 1.2,
    cs -> exclusive_jets(cs; njets = 4), "exclusive njets beta=1.2 gamma=1.2 R=0.8"
)
run_reco_test(valencia_exclusive4)

valencia_exclusive_d500 = ComparisonTestValencia(
    events_file_ee,
    joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-exclusive-d500-beta1.2-gamma1.2-R0.8-eeH.json.zst"),
    JetAlgorithm.Valencia, RecoStrategy.N2Plain, 1.2, 0.8, 1.2,
    cs -> exclusive_jets(cs; dcut = 500.0), "exclusive dcut=500 beta=1.2 gamma=1.2 R=0.8"
)
run_reco_test(valencia_exclusive_d500)

Reference outputs for R=1.0, beta=1.0, gamma=1.0
valencia_inclusive_b1g1 = ComparisonTestValencia(
    events_file_ee,
    joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-inclusive20-beta1.0-gamma1.0-R1.0-eeH.json.zst"),
    JetAlgorithm.Valencia, RecoStrategy.N2Plain, 1.0, 1.0, 1.0,
    cs -> inclusive_jets(cs; ptmin = 20.0), "inclusive beta=1.0 gamma=1.0 R=1.0"
)
run_reco_test(valencia_inclusive_b1g1)

valencia_exclusive4_b1g1 = ComparisonTestValencia(
    events_file_ee,
    joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-exclusive4-beta1.0-gamma1.0-R1.0-eeH.json.zst"),
    JetAlgorithm.Valencia, RecoStrategy.N2Plain, 1.0, 1.0, 1.0,
    cs -> exclusive_jets(cs; njets = 4), "exclusive njets beta=1.0 gamma=1.0 R=1.0"
)
run_reco_test(valencia_exclusive4_b1g1)

valencia_exclusive_d500_b1g1 = ComparisonTestValencia(
    events_file_ee,
    joinpath(@__DIR__, "data", "jet-collections-fastjet-valencia-exclusive-d500-beta1.0-gamma1.0-R1.0-eeH.json.zst"),
    JetAlgorithm.Valencia, RecoStrategy.N2Plain, 1.0, 1.0, 1.0,
    cs -> exclusive_jets(cs; dcut = 500.0), "exclusive dcut=500 beta=1.0 gamma=1.0 R=1.0"
)
run_reco_test(valencia_exclusive_d500_b1g1)

