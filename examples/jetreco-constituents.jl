# # Jet Reconstruction Constituents Example
#
# Perform a simple reconstruction example and show how to retrieve constituent jets.
#
# N.B. currently you must use the `jet-constituents` branch of `JetReconstruction`.
using JetReconstruction
using LorentzVectorHEP

input_file = joinpath(dirname(pathof(JetReconstruction)), "..", "test", "data",
                      "events.hepmc3.gz")
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no], p = 1, R = 1.0)

# Retrieve the exclusive pj_jets, but as `PseudoJet` types
pj_jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

# Get the constituents of the first jet
my_constituents = JetReconstruction.constituents(pj_jets[1], cluster_seq)

println("Constituents of jet number $(event_no):")
for c in my_constituents
    println(" $c")
end

# Now show how to convert to LorentzVectorCyl:
println("\nConstituents of jet number $(event_no) as LorentzVectorCyl:")
for c in my_constituents
    println(" $(LorentzVectorCyl(JetReconstruction.pt(c), JetReconstruction.rapidity(c), JetReconstruction.phi(c), JetReconstruction.mass(c)))")
end
