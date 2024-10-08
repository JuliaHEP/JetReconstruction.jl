# This example is suitable for running in the REPL
using EDM4hep
using EDM4hep.RootIO
using JetReconstruction

input_file = joinpath("/", "Users", "graemes", "code", "EDM4hepJets", "data",
                      "events_196755633.root")
reader = RootIO.Reader(input_file)
events = RootIO.get(reader, "events")

evt = events[1]

recps = RootIO.get(reader, evt, "ReconstructedParticles")

cs = jet_reconstruct(recps; algorithm = JetAlgorithm.Durham)
for jet in exclusive_jets(cs; njets = 2, T = EEjet)
    println(jet)
end
