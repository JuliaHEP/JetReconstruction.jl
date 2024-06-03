## Jet visualisation

module JetVisualisation

using JetReconstruction
using CairoMakie

println("JetVisualisation module has been loaded")

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
function JetReconstruction.jetsplot(objects, cs::ClusterSequence; barsize_phi = 0.1, barsize_eta = 0.1, colormap = :glasbey_hv_n256, Module = CairoMakie)
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
function JetReconstruction.jetsplot(objects, idx_arrays; barsize_phi = 0.1, barsize_eta = 0.1, colormap = :glasbey_hv_n256, Module = CairoMakie)
    cs = fill(0, length(objects)) # colours
    for i in 1:length(idx_arrays), j in idx_arrays[i]
        cs[j] = i
    end

    pts = sqrt.(JetReconstruction.pt2.(objects))

    Module.meshscatter(
        Module.Point3f.(JetReconstruction.phi.(objects), JetReconstruction.rapidity.(objects), 0pts);
        color = cs,
        markersize = Module.Vec3f.(barsize_phi, barsize_eta, pts),
        colormap = colormap,
        marker = Module.Rect3f(Module.Vec3f(0), Module.Vec3f(1)),
        figure = (size = (700, 600),),
        axis = (
            type = Module.Axis3, perspectiveness = 0.5, azimuth = 2.6, elevation = 0.5,
            xlabel = "ϕ", ylabel = "η", zlabel = "kt",
            limits = (nothing, nothing, nothing, nothing, 0, findmax(pts)[1] + 10),
        ),
        shading = NoShading,
    )
end

end
