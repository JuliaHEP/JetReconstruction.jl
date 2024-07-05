## Jet visualisation

module JetVisualisation

using JetReconstruction
using Makie

"""
    get_all_ancestors(idx, cs::ClusterSequence)
"""
function get_all_ancestors(idx, cs::ClusterSequence)
    if cs.history[idx].parent1 == JetReconstruction.NonexistentParent
        return [cs.history[idx].jetp_index]
    else
        branch1 = get_all_ancestors(cs.history[idx].parent1, cs)
        cs.history[idx].parent2 == JetReconstruction.BeamJet && return branch1
        branch2 = get_all_ancestors(cs.history[idx].parent2, cs)
        return [branch1; branch2]
    end
end

"""
    jetsplot(objects, cs::ClusterSequence; barsize_phi=0.1, barsize_eta=0.1, colormap=:glasbey_hv_n256, Module=Main)

Plots a 3d bar chart that represents jets. Takes `objects`, an array of objects to display (should be the same array you have passed to `jet_reconstruct` to get the `cs::ClusterSequence`), and the `cs::ClusterSequence` itself as arguments.

Optional arguments:
`barsize_phi::Real` — width of a bar along the ϕ axis;
`barsize_eta::Real` — width of a bar along the η axis;
`colormap::Symbol` — Makie colour map;
`Module` — the module where you have your Makie (see below);
```
# example
using CairoMakie # use any other Makie that you have here
jetsplot([object1, object2, object3], cluster_sequence_I_got_from_jet_reconstruct; Module=CairoMakie)
```

This function needs `Makie.jl` to work. You should install and import/use a specific backend yourself. `jetsplot` works with `CairoMakie`, `WGLMakie`, `GLMakie`, etc. Additionally, you can specify the module where you have your `Makie` explicitly:
```
import CairoMakie
jetsplot(my_objects, cs, Module=CairoMakie)

import GLMakie
jetsplot(my_objects, cs, Module=GLMakie)

using WGLMakie
jetsplot(my_objects, cs, Module=Main) #default
```
"""
function JetReconstruction.jetsplot(objects, cs::ClusterSequence; barsize_phi = 0.1,
                                    barsize_eta = 0.1, colormap = :glasbey_hv_n256,
                                    Module = Makie)
    idx_arrays = Vector{Int}[]
    for elt in cs.history
        elt.parent2 == JetReconstruction.BeamJet || continue
        push!(idx_arrays, get_all_ancestors(elt.parent1, cs))
    end

    jetsplot(objects, idx_arrays; barsize_phi, barsize_eta, colormap, Module)
end

"""
`jetsplot(objects, idx_arrays; barsize_phi=0.1, barsize_eta=0.1, colormap=:glasbey_hv_n256, Module=Main)`

Plots a 3d bar chart that represents jets. Takes an `objects` array of objects to display and `idx_arrays`, an array of arrays with indeces, where `idx_arrays[i]` gives indeces of `objects` that form the jet number `i`. This function's signature might not be the most practical for the current version of the JetReconstruction.jl package, as it has been written during the early stage of development. There is now an overload of it that takes a `ClusterSequence` object as its argument.

Optional arguments:
`barsize_phi::Real` — width of a bar along the ϕ axis;
`barsize_eta::Real` — width of a bar along the η axis;
`colormap::Symbol` — Makie colour map;
`Module` — the module where you have your Makie (see below);
```
# example
using CairoMakie # use any other Makie that you have here

jetsplot([object1, object2, object3], [[1], [2, 3]])
```
The example above plots `object1` as a separate jet in one colour and `object2` and `object3` together in another colour.

This function needs `Makie.jl` to work. You should install and import/use a specific backend yourself. `jetsplot` works with `CairoMakie`, `WGLMakie`, `GLMakie`, etc. Additionally, you can specify the module where you have your `Makie` explicitly:
```
import CairoMakie
jetsplot(my_objects, my_colour_arrays, Module=CairoMakie)

import GLMakie
jetsplot(my_objects, my_colour_arrays, Module=GLMakie)

using WGLMakie
jetsplot(my_objects, my_colour_arrays, Module=Main) #default
```
"""
function JetReconstruction.jetsplot(objects, idx_arrays; barsize_phi = 0.1,
                                    barsize_eta = 0.1, colormap = :glasbey_hv_n256,
                                    Module = Main)
    cs = fill(0, length(objects)) # colours
    for i in 1:length(idx_arrays), j in idx_arrays[i]
        cs[j] = i
    end

    pts = sqrt.(JetReconstruction.pt2.(objects))

    Module.meshscatter(Module.Point3f.(JetReconstruction.phi.(objects),
                                       JetReconstruction.rapidity.(objects), 0pts);
                       color = cs,
                       markersize = Module.Vec3f.(barsize_phi, barsize_eta, pts),
                       colormap = colormap,
                       marker = Module.Rect3f(Module.Vec3f(0), Module.Vec3f(1)),
                       figure = (size = (700, 600),),
                       axis = (type = Module.Axis3, perspectiveness = 0.5, azimuth = 2.6,
                               elevation = 0.5,
                               xlabel = "ϕ", ylabel = "η", zlabel = "kt",
                               limits = (nothing, nothing, nothing, nothing, 0,
                                         findmax(pts)[1] + 10)),
                       shading = NoShading)
end

function JetReconstruction.jetsplot(cs::ClusterSequence,
                                    reco_state::Dict{Int,
                                                     JetReconstruction.JetWithAncestors};
                                    barsize_phi = 0.1,
                                    barsize_y = 0.1, colormap = :glasbey_category10_n256,
                                    Module = Makie)
    # Setup the marker as a square object
    jet_plot_marker = Rect3f(Vec3f(0), Vec3f(1))

    # Get the jet variables to plot
    phis = [JetReconstruction.phi(x.self) for x in values(reco_state)]
    ys = [JetReconstruction.rapidity(x.self) for x in values(reco_state)]
    pts = [JetReconstruction.pt(x.self) for x in values(reco_state)]

    # The core points to plot are on the pt=0 axis, with the marker size
    # scaled up to the p_T of the jet
    jet_plot_points = Point3f.(phis .- (barsize_phi / 2), ys .- (barsize_y / 2), 0pts)
    jet_plot_marker_size = Vec3f.(barsize_phi, barsize_y, pts)

    # Colours are defined from the rank of the ancestor jets
    jet_plot_colours = [x.jet_rank for x in values(reco_state)]

    # Limits for rapidity and p_T are set to the maximum values in the data
    # (For ϕ the limits are (0, 2π))
    min_rap = max_rap = max_pt = 0.0
    for jet in cs.jets
        min_rap = min(min_rap, JetReconstruction.rapidity(jet))
        max_rap = max(max_rap, JetReconstruction.rapidity(jet))
        max_pt = max(max_pt, JetReconstruction.pt(jet))
    end

    fig, ax, plt_obj = Module.meshscatter(jet_plot_points;
                                          markersize = jet_plot_marker_size,
                                          marker = jet_plot_marker,
                                          colormap = colormap,
                                          color = jet_plot_colours,
                                          colorrange = (1, 256),
                                          figure = (size = (700, 600),),
                                          axis = (type = Axis3, perspectiveness = 0.5,
                                                  azimuth = 2.7,
                                                  elevation = 0.5,
                                                  xlabel = L"\phi", ylabel = L"y",
                                                  zlabel = L"p_T",
                                                  limits = (0, 2π, min_rap - 0.5,
                                                            max_rap + 0.5, 0, max_pt + 10)),
                                          shading = NoShading)
    fig, ax, plt_obj
end

"""
    animatereco(cs::ClusterSequence, filename;
                barsize_phi = 0.1,
                barsize_y = 0.1,
                colormap = :glasbey_category10_n256,
                perspective = 0.5,
                azimuth = 2.7,
                elevation = 0.5,
                framerate = 5,
                ancestors = false,
                Module = Makie)

Animate the jet reconstruction process and save it as a video file.

# Arguments
- `cs::ClusterSequence`: The cluster sequence object containing the jets.
- `filename`: The name of the output video file.

# Optional Arguments
- `barsize_phi=0.1`: The size of the bars in the phi direction.
- `barsize_y=0.1`: The size of the bars in the y direction.
- `colormap=:glasbey_category10_n256`: The colormap to use for coloring the
  jets.
- `perspective=0.5`: The perspective of the plot.
- `azimuth=2.7`: The azimuth angle of the plot.
- `elevation=0.5`: The elevation angle of the plot.
- `framerate=5`: The framerate of the output video.
- `ancestors=false`: Whether to include ancestors of the jets in the animation.
  When `true` the ancestors of the jets will be plotted as well, as height zero
  bars, with the same colour as the jet they are ancestors of.
- `Module`: The plotting module to use. Default is `Makie`.

For `perspective`, `azimuth`, and `elevation`, a single value can be passed for
a fixed viewpoint, or a tuple of two values for a changing viewpoint. The
viewpoint will then change linearly between the two values over the course of the
animation.

# Returns
- `fig`: The figure object representing the final fram.

"""
function JetReconstruction.animatereco(cs::ClusterSequence, filename;
                                       barsize_phi = 0.1,
                                       barsize_y = 0.1,
                                       colormap = :glasbey_category10_n256,
                                       perspective::Union{Real, Tuple{Real, Real}} = 0.5,
                                       azimuth::Union{Real, Tuple{Real, Real}} = 2.7,
                                       elevation::Union{Real, Tuple{Real, Real}} = 0.5,
                                       framerate = 5, ancestors = false,
                                       Module = Makie)
    # Setup the marker as a square object
    jet_plot_marker = Rect3f(Vec3f(0), Vec3f(1))

    # Get the number of meaningful reconstruction steps and rank the initial
    # particles by p_T
    merge_steps = JetReconstruction.merge_steps(cs)
    jet_ranks = JetReconstruction.jet_ranks(cs)

    # End point of the catagorical color map
    # (How to get this programmatically from the CM Symbol?)
    colormap_end = 256

    # Get the reconstruction state at each meaningful iteration
    # And calculate the plot parameters
    all_jet_plot_points = Vector{Vector{Point3f}}()
    all_jet_plot_marker_size = Vector{Vector{Vec3f}}()
    all_jet_plot_colours = Vector{Vector{Int}}()
    for step in 0:merge_steps
        reco_state = JetReconstruction.reco_state(cs, jet_ranks; iteration = step)
        phis = [JetReconstruction.phi(x.self) for x in values(reco_state)]
        ys = [JetReconstruction.rapidity(x.self) for x in values(reco_state)]
        pts = [JetReconstruction.pt(x.self) for x in values(reco_state)]
        push!(all_jet_plot_points,
              Point3f.(phis .- (barsize_phi / 2), ys .- (barsize_y / 2), 0pts))
        push!(all_jet_plot_marker_size, Vec3f.(barsize_phi, barsize_y, pts))
        push!(all_jet_plot_colours, [mod1(x.jet_rank, colormap_end) for x in values(reco_state)])
        if ancestors
            for jet_entry in values(reco_state)
                for ancestor in jet_entry.ancestors
                    ancestor_jet = cs.jets[ancestor]
                    push!(all_jet_plot_points[end],
                          Point3f(JetReconstruction.phi(ancestor_jet) - (barsize_phi / 2),
                                  JetReconstruction.rapidity(ancestor_jet) -
                                  (barsize_y / 2),
                                  0.0))
                    push!(all_jet_plot_marker_size[end],
                          Vec3f(barsize_phi, barsize_y, 0.001))
                    push!(all_jet_plot_colours[end], mod1(jet_entry.jet_rank, colormap_end))
                end
            end
        end
    end

    # Keep plot limits constant
    min_rap = max_rap = max_pt = 0.0
    for jet in cs.jets
        min_rap = min(min_rap, JetReconstruction.rapidity(jet))
        max_rap = max(max_rap, JetReconstruction.rapidity(jet))
        max_pt = max(max_pt, JetReconstruction.pt(jet))
    end

    # Setup an observable for the iteration number
    it_obs = Observable(0)
    jet_plot_points_obs = @lift all_jet_plot_points[$it_obs + 1]
    jet_plot_marker_size_obs = @lift all_jet_plot_marker_size[$it_obs + 1]
    jet_plot_colours_obs = @lift all_jet_plot_colours[$it_obs + 1]

    # We may want to have a shifting viewpoint
    azimuth_axis = typeof(azimuth) <: Tuple ?
                   @lift(azimuth[1]+$it_obs / merge_steps * (azimuth[2] - azimuth[1])) :
                   azimuth
    elevation_axis = typeof(elevation) <: Tuple ?
                     @lift(elevation[1]+$it_obs / merge_steps *
                                        (elevation[2] - elevation[1])) : elevation
    perspective_axis = typeof(perspective) <: Tuple ?
                       @lift(perspective[1]+$it_obs / merge_steps *
                                            (perspective[2] - perspective[1])) : perspective

    ax = (type = Axis3,
          xlabel = L"\phi", ylabel = L"y", zlabel = L"p_T",
          limits = (0, 2π, min_rap - 0.5, max_rap + 0.5, 0, max_pt + 10),
          perspectiveness = perspective_axis, azimuth = azimuth_axis,
          elevation = elevation_axis)
    fig = Module.meshscatter(jet_plot_points_obs;
                             markersize = jet_plot_marker_size_obs,
                             marker = jet_plot_marker,
                             colormap = colormap,
                             color = jet_plot_colours_obs,
                             colorrange = (1, colormap_end),
                             figure = (size = (800, 600),),
                             axis = ax,
                             shading = NoShading)
    record(fig, filename, 0:merge_steps; framerate = framerate) do iteration
        it_obs[] = iteration
    end
    fig
end

end
