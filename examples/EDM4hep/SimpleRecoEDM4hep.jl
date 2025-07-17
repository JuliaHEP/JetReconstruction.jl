# This example is suitable for running in the REPL
using EDM4hep
using EDM4hep.RootIO
using JetReconstruction

input_file = joinpath("/", "Users", "graemes", "code", "EDM4hepJets", "data",
                      "events_196755633.root")
reader = RootIO.Reader(input_file)
events = RootIO.get(reader, "events");

evt = events[1];

recps = RootIO.get(reader, evt, "ReconstructedParticles")

# Reconstruct and print the jets
cs = jet_reconstruct(recps; algorithm = JetAlgorithm.Durham)
dijets = exclusive_jets(cs, EEJet; njets = 2)
for jet in dijets
    println(jet)
end

# Get constituents
for (i, jet) in enumerate(dijets)
    my_constituent_indexes = constituent_indexes(jet, cs)
    println("Jet $i constituents: $my_constituent_indexes")
end
