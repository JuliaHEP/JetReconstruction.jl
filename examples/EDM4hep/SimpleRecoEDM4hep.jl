# This example is suitable for running in the REPL
using EDM4hep
using EDM4hep.RootIO
using JetReconstruction

input_file = joinpath(pkgdir(EDM4hep), "examples", "ttbar_edm4hep_digi.root")
reader = RootIO.Reader(input_file)

events = RootIO.get(reader, "events");

evt = events[1];

recps = RootIO.get(reader, evt, "MCParticle")

# Reconstruct and print the jets
cs = jet_reconstruct(recps; algorithm = JetAlgorithm.AntiKt)
jets = inclusive_jets(cs, PseudoJet; ptmin = 5.0)
for jet in jets
    println(jet)
end

# Get constituents
for (i, jet) in enumerate(jets)
    my_constituent_indexes = constituent_indexes(jet, cs)
    println("Jet $i constituents: $my_constituent_indexes")
end
