using JetReconstruction
using CairoMakie
using StatsBase

"""
generate_average_lund_image(njets::Int, delta_array::Vector{Vector{Real}}, 
                            kt_array::Vector{Vector{Real}}; xrange::Tuple{Real}, 
                            yrange::Tuple{Real}, bins::Int) -> (xgrid, ygrid, avg_image)

Computes an average Lund image from a set of jets by binning the (log(1/ΔR), log(kt))
coordinates into a fixed-size 2D histogram. Each jet's histogram is normalized and
then averaged across all jets.

Arguments:
- `njets`: Number of jets
- `delta_array`: Vector of vectors, where each inner vector contains ΔR values for a jet
- `kt_array`: Vector of vectors, where each inner vector contains kt values for a jet
- `xrange`: Tuple defining the x-axis (log(1/ΔR)) range
- `yrange`: Tuple defining the y-axis (log(kt)) range
- `bins`: Number of bins along each axis

Returns:
- A tuple `(xgrid, ygrid, avg_image)` where `xgrid` and `ygrid` are coordinate axes labels and
  `avg_image` is the averaged 2D histogram as a matrix
"""
function generate_average_lund_image(njets::Int, delta_array::Vector{Vector{Real}},
                                     kt_array::Vector{Vector{Real}}; xrange = (0.0, 4.0),
                                     yrange = (-5.0, 7.0), bins = 25)
    xmin, xmax = xrange
    ymin, ymax = yrange

    x_width = (xmax - xmin) / bins
    y_width = (ymax - ymin) / bins

    total_res = []

    for i in 1:njets
        Xind = ceil.((delta_array[i] .- xmin) ./ x_width)
        Yind = ceil.((kt_array[i] .- ymin) ./ y_width)

        res = zeros(Float32, bins, bins)
        L1norm = 0.0

        for i in 1:length(Xind)
            x = Int(Xind[i])
            y = Int(Yind[i])
            if (maximum([x, y]) < bins && minimum([x, y]) >= 1)
                res[x, y] += 1
                L1norm += 1
            end
        end

        if L1norm > 0
            res[1, :] .= res[1, :] ./ L1norm
            push!(total_res, res)
        end
    end

    avg_res = mean(total_res, dims = 1)[1]
    x = range(xmin, xmax; length = bins)
    y = range(ymin, ymax; length = bins)

    return (x, y, avg_res)
end

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
    jets = sort!(inclusive_jets(cluster_seq, PseudoJet; ptmin = 10.0),
                 by = JetReconstruction.pt2, rev = true)

    @info "Generating Primary Lund Emissions for $(length(jets)) jets for Event $(event_no):"
    for (ijet, jet) in enumerate(jets)
        lundvars = generate_lund_emissions(jet, cluster_seq)
        println("- Jet $(ijet) has $(length(lundvars)) emissions in lund plane")
        kt = lund_kt.(lundvars)
        Δ = lund_delta.(lundvars)

        push!(lundX, -log.(Δ))
        push!(lundY, log.(kt))
    end
end

njets = length(lundX)
x, y, avg_res = generate_average_lund_image(njets, lundX, lundY)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "ln(R/Δ)", ylabel = "ln(kt/GeV)",
          title = "Average Lund Image")
hm = heatmap!(ax, x, y, avg_res; colormap = :viridis, colorrange = extrema(avg_res))
Colorbar(fig[1, 2], hm; label = "")

display(fig)

save("lund-gen.png", fig)
