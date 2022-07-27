## Jet visualisation
# not a submodule

"""
`jetsplot(objects, idx_arrays; barsize=0.1, colormap=:glasbey_hv_n256, Module=Main)`

Plots a 3d bar chart that represents jets. Takes an `objects` array of objects to display and `idx_arrays`, an array of arrays with indeces, where `idx_arrays[i]` gives indeces of `objects` that form the jet number `i`.

Optional arguments:
`barsize::Real` — the width of the individual bars;
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
function jetsplot(objects, idx_arrays; barsize=0.1, colormap=:glasbey_hv_n256, Module=Main)
	cs = fill(0, length(objects))
	for i in 1:length(idx_arrays), j in idx_arrays[i]
		cs[j] = i
	end

	pts = pt.(objects)

	Module.meshscatter(
		Module.Point3f.(phi.(objects), eta.(objects), 0pts);
	  	color = cs,
		markersize = Module.Vec3f.(barsize, barsize, pts),
		colormap = colormap,
		marker = Module.Rect3f(Module.Vec3f(0), Module.Vec3f(1)),
	 	figure = (resolution=(700,600),),
		axis = (
			type = Module.Axis3, perspectiveness = 0.5, azimuth = 2.6, elevation=0.5,
                        xlabel = "ϕ", ylabel = "η", zlabel = "kt",
		        limits = (nothing, nothing, nothing, nothing, 0, findmax(pts)[1]+10)
		),
	    shading=false
	)
end
