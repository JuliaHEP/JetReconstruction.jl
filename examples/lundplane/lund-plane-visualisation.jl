using JetReconstruction
using CairoMakie
using StatsBase

input_file = joinpath(dirname(pathof(JetReconstruction)),
                      "..", "test", "data", "events.pp13TeV.hepmc3.zst")
events = read_final_state_particles(input_file)

N = length(events)
lundX = Vector{Vector{Real}}()
lundY = Vector{Vector{Real}}()

lund_kt(p) = p.kt
lund_delta(p) = p.delta

# Event to pick
for event_no in 1:N
    cluster_seq = jet_reconstruct(events[event_no]; algorithm = JetAlgorithm.AntiKt,
                                  R = 1.0)
    jets = sort!(inclusive_jets(cluster_seq; T = PseudoJet, ptmin = 10.0),
                 by = JetReconstruction.pt2, rev = true)

    @info "Generating Lund Jets for $(length(jets)) jets for Event $(event_no):"
    for (ijet, jet) in enumerate(jets)
        lundvars = generate_lund_projection(jet, cluster_seq)
        println("- Jet $(ijet) has $(length(lundvars)) projections in lund plane")
        kt = lund_kt.(lundvars)
        Δ = lund_delta.(lundvars)

        push!(lundX, -log.(Δ))
        push!(lundY, log.(kt))
    end
end

njets = length(lundX)
println(typeof(lundX))
x, y, avg_res = generate_average_lund_image(njets, lundX, lundY)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "ln(R/Δ)", ylabel = "ln(kt/GeV)",
          title = "Average Lund Image")
hm = heatmap!(ax, x, y, avg_res; colormap = :viridis, colorrange = extrema(avg_res))
Colorbar(fig[1, 2], hm; label = "")

display(fig)

save("lund-gen.png", fig)
