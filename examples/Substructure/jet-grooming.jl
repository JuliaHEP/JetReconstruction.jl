#! /usr/bin/env julia
using JetReconstruction

input_file = joinpath(dirname(pathof(JetReconstruction)),
                      "..", "test", "data", "events.pp13TeV.hepmc3.gz")
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no], p = 0, R = 1.0)
jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

r = 0.3     # recluster radius
n = 3       # number of hard jets to consider

filter = JetFilter(r, n)

@info "Jet Filtering: recluster radius = $r, hard subjets to consider = $n"
for jet in jets
    filtered = jet_filtering(jet, cluster_seq, filter)

    println("Original jet: pt = $(JetReconstruction.pt(jet)), rap = $(JetReconstruction.rapidity(jet)), phi = $(JetReconstruction.phi(jet)), E = $(jet.E)")
    println("Filtered jet: pt = $(JetReconstruction.pt(filtered)), rap = $(JetReconstruction.rapidity(filtered)), phi = $(JetReconstruction.phi(filtered)), E = $(filtered.E)\n")
end

r = 0.3                     # recluster radius
f = 0.3                     # trim fraction
m = JetAlgorithm.CA         # recluster method

trim = JetTrim(r, f, m)

@info "Jet Trimming: recluster radius = $r, trim fraction = $f, recluster method = $m"
for jet in jets
    trimmed = jet_trimming(jet, cluster_seq, trim)

    println("Original jet: pt = $(JetReconstruction.pt(jet)), rap = $(JetReconstruction.rapidity(jet)), phi = $(JetReconstruction.phi(jet)), E = $(jet.E)")
    println("Trimmed jet: pt = $(JetReconstruction.pt(trimmed)), rap = $(JetReconstruction.rapidity(trimmed)), phi = $(JetReconstruction.phi(trimmed)), E = $(trimmed.E)\n")
end
