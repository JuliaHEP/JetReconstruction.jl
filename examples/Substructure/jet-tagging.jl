#! /usr/bin/env julia
using JetReconstruction

input_file = joinpath(dirname(pathof(JetReconstruction)),
                      "..", "test", "data", "events.pp13TeV.hepmc3.gz")
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no], p = 0, R = 1.0)
jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

μ = 0.67        # jet mass ratio
y = 0.09        # symmetry cut

MDtagger = MassDropTagger(μ, y)

@info "Mass Drop Tagging: μ = $μ, y = $y"
for jet in jets
    tagged = mass_drop(jet, cluster_seq, MDtagger)
    println("Original jet: pt = $(JetReconstruction.pt(jet)), rap = $(JetReconstruction.rapidity(jet)), phi = $(JetReconstruction.phi(jet)), E = $(jet.E)")
    println("Tagged jet: pt = $(JetReconstruction.pt(tagged)), rap = $(JetReconstruction.rapidity(tagged)), phi = $(JetReconstruction.phi(tagged)), E = $(tagged.E)\n")
end

z = 0.1         # soft drop threshold
b = 2.0         # angular exponent

SDtagger = SoftDropTagger(z, b)

@info "Soft Drop Tagging: recluster radius = $(SDtagger.cluster_rad), zcut = $z, b = $b"
for jet in jets
    tagged = soft_drop(jet, cluster_seq, SDtagger)
    println("Original jet: pt = $(JetReconstruction.pt(jet)), rap = $(JetReconstruction.rapidity(jet)), phi = $(JetReconstruction.phi(jet)), E = $(jet.E)")
    println("Tagged jet: pt = $(JetReconstruction.pt(tagged)), rap = $(JetReconstruction.rapidity(tagged)), phi = $(JetReconstruction.phi(tagged)), E = $(tagged.E)\n")
end
