using JetReconstruction

input_file = joinpath(dirname(pathof(JetReconstruction)),
                      "..", "test", "data", "events.pp13TeV.hepmc3.zst")
# input_file = "/Users/sattwamoghosh/IISER/Sem7/PH4112/hepmc_file/pp_to_tt.hepmc"
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no]; algorithm = JetAlgorithm.AntiKt, R = 1.0)
jets = sort!(inclusive_jets(cluster_seq; T = PseudoJet, ptmin = 10.0),
             by = JetReconstruction.pt2, rev = true)

@info "Generating Lund Jets for $(length(jets)) jets for Event $(event_no):"
for (ijet, jet) in enumerate(jets)
    lundvars = generate_lund_projection(jet, cluster_seq)
    println("- Jet $(ijet) has $(length(lundvars)) projections in lund plane")
end
