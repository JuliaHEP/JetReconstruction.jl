## Jet visualisation

using PyPlot

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
    phidata = phi.(objects)
    etadata = eta.(objects)
    ktdata = pt.(objects)

    fig = figure()
    ax = fig.add_subplot(projection="3d")

    ax.set_ylabel("η")
    ax.set_xlabel("ϕ")
    ax.set_zlabel("kt")
    for idx in idx_arrays
        ax.bar3d(phidata[idx], etadata[idx], 0ktdata[idx], barsize, barsize, ktdata[idx])
    end

    PyPlot.gcf()
end
