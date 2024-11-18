# Dummy code to call and cache compilation of all C-bindings
# The numeric values doesn't make sense and are used only to select correct dispatch

using JetReconstruction
using JetReconstruction.C_JetReconstruction: jetreconstruction_PseudoJet_init,
                                             C_ClusterSequence,
                                             jetreconstruction_ClusterSequence_free_members_,
                                             C_JetsResult,
                                             jetreconstruction_JetsResult_free_members_,
                                             jetreconstruction_jet_reconstruct,
                                             jetreconstruction_inclusive_jets,
                                             jetreconstruction_exclusive_jets_njets,
                                             jetreconstruction_exclusive_jets_dcut

pseudoJets_len = Csize_t(2)
pseudoJets_ptr = Ptr{PseudoJet}(Libc.malloc(pseudoJets_len * sizeof(PseudoJet)))

jetreconstruction_PseudoJet_init(pseudoJets_ptr, 0.0, 1.0, 2.0, 3.0)
jetreconstruction_PseudoJet_init(pseudoJets_ptr + sizeof(PseudoJet), 1.0, 2.0, 3.0, 4.0)

strategy = RecoStrategy.Best
algorithm = JetAlgorithm.CA
R = 2.0

clustersequence_ptr = Ptr{C_ClusterSequence{PseudoJet}}(Libc.malloc(sizeof(C_ClusterSequence{PseudoJet})))
jetreconstruction_jet_reconstruct(pseudoJets_ptr, pseudoJets_len, algorithm, R,
                                  RecoStrategy.Best,
                                  clustersequence_ptr)

results = C_JetsResult{PseudoJet}(C_NULL, 0)
results_ptr = Base.unsafe_convert(Ptr{C_JetsResult{PseudoJet}}, Ref(results))
jetreconstruction_exclusive_jets_njets(clustersequence_ptr, Csize_t(2), results_ptr)
jetreconstruction_JetsResult_free_members_(results_ptr)

jetreconstruction_exclusive_jets_dcut(clustersequence_ptr, 1.0, results_ptr)
jetreconstruction_JetsResult_free_members_(results_ptr)

jetreconstruction_inclusive_jets(clustersequence_ptr, 0.0, results_ptr)
jetreconstruction_JetsResult_free_members_(results_ptr)

jetreconstruction_ClusterSequence_free_members_(clustersequence_ptr)
Libc.free(pseudoJets_ptr)
Libc.free(clustersequence_ptr)
