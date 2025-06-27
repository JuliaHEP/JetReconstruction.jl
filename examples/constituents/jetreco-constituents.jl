# # Jet Reconstruction Constituents Example
#
# Perform a simple reconstruction example and show how to retrieve constituent jets.
using JetReconstruction
using LorentzVectorHEP
using Logging

logger = ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

input_file = joinpath(dirname(pathof(JetReconstruction)), "..", "test", "data",
                      "events.pp13TeV.hepmc3.zst")
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no]; algorithm = JetAlgorithm.Kt, R = 1.0)

# Retrieve the exclusive pj_jets, but as `PseudoJet` types
pj_jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

# Get the constituents of the first jet
my_constituents = JetReconstruction.constituents(pj_jets[1], cluster_seq)

println("Constituents of jet number $(event_no):")
for c in my_constituents
    println(" $c")
end

# Just retrieve the indexes of the constituents
my_constituent_indexes = constituent_indexes(pj_jets[1], cluster_seq)
println("\nConstituent indexes for jet number $(event_no): $my_constituent_indexes")
for i in my_constituent_indexes
    println("  Constituent jet $i: $(events[1][i])")
end
