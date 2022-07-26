## Jet visualisation

using CairoMakie
#CairoMakie.activate!(type = "svg")

"""
`jetsplot(objects, idx_arrays; barsize=0.1)`

Plots a 3d bar chart that represents jets. Takes an `objects` array of objects to display and `idx_arrays`, an array of arrays with indeces, where `idx_arrays[i]` gives indeces of objects that form the jet number `i`.
```
# example
jetsplot([object1, object2, object3], [[1], [2, 3]])
```
The example above plots `object1` as a separate jet in one colour and `object2` and `object3` together in another colour.
"""
function jetsplot(objects, idx_arrays; barsize=0.1)
	cs = fill(0, length(objects))
	for i in 1:length(idx_arrays), j in idx_arrays[i]
		cs[j] = i
	end

	pts = pt.(objects)

	meshscatter(
		Point3f.(phi.(objects), eta.(objects), 0pts);
	  	color = cs,
		markersize = Vec3f.(barsize, barsize, pts),
		colormap = :Spectral_7,
		marker = Rect3f(Vec3f(0), Vec3f(1)),
	 	figure = (resolution=(700,600),),
		axis = (
			type = Axis3, perspectiveness = 0.5, azimuth = 2.6, elevation=0.5,
                        xlabel = "ϕ", ylabel = "η", zlabel = "kt",
		        limits = (nothing, nothing, nothing, nothing, 0, findmax(pts)[1]+10)
		),
	    shading=false
	)
end
