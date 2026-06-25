# Dummy code to call and cache compilation of all C-bindings
# The numeric values don't make sense and are used only to select correct dispatch

using JetReconstruction: PseudoJet, RecoStrategy, JetAlgorithm, RecombinationScheme
using JetReconstruction.C_JetReconstruction: jetreconstruction_PseudoJet_init,
                                             C_ClusterSequence,
                                             jetreconstruction_ClusterSequence_free_members_,
                                             C_JetsResult,
                                             jetreconstruction_JetsResult_free_members_,
                                             jetreconstruction_jet_reconstruct,
                                             jetreconstruction_inclusive_jets,
                                             jetreconstruction_exclusive_jets_njets,
                                             jetreconstruction_exclusive_jets_dcut,
                                             StatusCode

function assert_ok(ret_val)
    @assert StatusCode.T(ret_val) == StatusCode.OK
end

R = 2.0
jets_len = Csize_t(2)
power = 0.5

### pp jets ###

jets_ptr = Ptr{PseudoJet}(Libc.malloc(jets_len * sizeof(PseudoJet)))
@assert jets_ptr != C_NULL
clustersequence_ptr = Ptr{C_ClusterSequence{PseudoJet}}(Libc.malloc(sizeof(C_ClusterSequence{PseudoJet})))
@assert clustersequence_ptr != C_NULL
results_ptr = Ptr{C_JetsResult{PseudoJet}}(Libc.malloc(sizeof(C_JetsResult{PseudoJet})))
@assert results_ptr != C_NULL

assert_ok(jetreconstruction_PseudoJet_init(jets_ptr, 0.0, 1.0, 2.0, 3.0, 1))
assert_ok(jetreconstruction_PseudoJet_init(jets_ptr + sizeof(PseudoJet), 1.0, 2.0, 3.0,
                                           4.0, 2))

for algorithm in [JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.GenKt],
    strategy in [RecoStrategy.N2Plain, RecoStrategy.N2Tiled, RecoStrategy.Best],
    recombination in [
        RecombinationScheme.EScheme,
        RecombinationScheme.PtScheme,
        RecombinationScheme.Pt2Scheme
    ]

    assert_ok(jetreconstruction_jet_reconstruct(jets_ptr, jets_len, algorithm, power, R,
                                                strategy, recombination,
                                                clustersequence_ptr))

    assert_ok(jetreconstruction_inclusive_jets(clustersequence_ptr, 0.0,
                                               results_ptr))
    jetreconstruction_JetsResult_free_members_(results_ptr)

    if algorithm != JetAlgorithm.AntiKt
        assert_ok(jetreconstruction_exclusive_jets_njets(clustersequence_ptr,
                                                         Csize_t(2),
                                                         results_ptr))
        jetreconstruction_JetsResult_free_members_(results_ptr)

        assert_ok(jetreconstruction_exclusive_jets_dcut(clustersequence_ptr, 1.0,
                                                        results_ptr))
        jetreconstruction_JetsResult_free_members_(results_ptr)

        jetreconstruction_ClusterSequence_free_members_(clustersequence_ptr)
    end
end
Libc.free(jets_ptr)
Libc.free(clustersequence_ptr)
Libc.free(results_ptr)
