#!/usr/bin/env julia

"""
Test script for Valencia algorithm implementation
Uses the first event from test/data/events.eeH.hepmc3.zst
Expected results from fastjet::contrib::ValenciaPlugin(1.2, 0.8)
"""

using JetReconstruction
using Printf

# Import the jet functions
import JetReconstruction: pt, rapidity, inclusive_jets, exclusive_jets, px, py, pz, energy

# Load the test event
eventfile = "test/data/events.eeH.hepmc3.zst"
events = read_final_state_particles(eventfile)

# Use the first event
event = events[1]

println("Testing Valencia algorithm with β=1.2, γ=0.8")
println("Input particles: $(length(event))")

# Run Valencia algorithm
β = 1.2
γ = 0.8
clusterseq = ee_genkt_algorithm(event, algorithm=JetAlgorithm.Valencia, p=β, γ=γ)

# Get inclusive jets (pT > 20 GeV)
inclusive_jets_result = inclusive_jets(clusterseq, ptmin=20.0)
println("\nHard jets (pT > 20 GeV) in inclusive clustering:")
println("Number of jets: $(length(inclusive_jets_result))")
for (i, jet) in enumerate(inclusive_jets_result)
    # Convert to PseudoJet for pt and rapidity functions
    pjet = PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
    @printf("pt = %.6g, rap = %.7g\n", pt(pjet), rapidity(pjet))
end

# Get exclusive jets (N=4)
exclusive_jets_n4 = exclusive_jets(clusterseq, njets=4)
println("\nHard jets in exclusive N=4 clustering:")
println("Number of jets: $(length(exclusive_jets_n4))")
for (i, jet) in enumerate(exclusive_jets_n4)
    # Convert to PseudoJet for pt and rapidity functions
    pjet = PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
    @printf("pt = %.6g, rap = %.7g\n", pt(pjet), rapidity(pjet))
end

# Get exclusive jets up to d=500
exclusive_jets_d500 = exclusive_jets(clusterseq, dcut=500.0)
println("\nHard jets in exclusive clustering up to d = 500:")
println("Number of jets: $(length(exclusive_jets_d500))")
for (i, jet) in enumerate(exclusive_jets_d500)
    # Convert to PseudoJet for pt and rapidity functions
    pjet = PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
    @printf("pt = %.6g, rap = %.7g\n", pt(pjet), rapidity(pjet))
end

println("\nExpected results from fastjet::contrib::ValenciaPlugin(1.2, 0.8):")
println("hard jets (pT > 20 GeV) in inclusive clustering")
println("pt = 122.875, rap = 0.0201274")
println("pt = 122.944, rap = -0.0140412")
println()
println("hard jets in exclusive N=4 clustering")
println("pt = 25.5392, rap = -0.0655267")
println("pt = 97.3956, rap = 0.0428278")
println("pt = 11.2882, rap = -0.122071")
println("pt = 112.065, rap = -0.0029899")
println()
println("hard jets in exclusive clustering up to d = 500")
println("pt = 122.944, rap = -0.0140412")
println("pt = 122.875, rap = 0.0201274")
