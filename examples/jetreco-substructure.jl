#! /usr/bin/env julia
using JetReconstruction

input_file = joinpath(dirname(pathof(JetReconstruction)), 
                            "..", "test", "data", "events.pp13TeV.hepmc3.gz")
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no], p = 0, R = 1.0)
jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

μ = 0.67
y = 0.09

tagger = MassDropTagger(μ, y)

@info "Mass Drop Tagging: μ = $μ, y = $y"
for jet in jets
    tagged = massDrop(jet, cluster_seq, tagger)
    println("Original jet: pt = $(JetReconstruction.pt(jet)), eta = $(JetReconstruction.eta(jet)), phi = $(JetReconstruction.phi(jet)), E = $(jet.E)")
    println("Tagged jet: pt = $(JetReconstruction.pt(tagged)), eta = $(JetReconstruction.eta(tagged)), phi = $(JetReconstruction.phi(tagged)), E = $(tagged.E)\n")
end

r = 0.3     # recluster radius
f = 0.3     # trim fraction
m = 0       # recluster method (0 = C/A, 1 = kT)

trim = Trim(r, f, m)

@info "Jet Trimming: recluster radius = $r, trim fraction = $f, recluster method = $(m == 0 ? "C/A" : "kT")"
for jet in jets
    trimmed = jetTrim(jet, cluster_seq, trim)

    println("Original jet: pt = $(JetReconstruction.pt(jet)), eta = $(JetReconstruction.eta(jet)), phi = $(JetReconstruction.phi(jet)), E = $(jet.E)")
    println("Trimmed jet: pt = $(JetReconstruction.pt(trimmed)), eta = $(JetReconstruction.eta(trimmed)), phi = $(JetReconstruction.phi(trimmed)), E = $(trimmed.E)\n")

end