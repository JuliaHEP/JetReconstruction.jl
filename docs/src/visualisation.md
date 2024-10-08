# Jet Visualisation Documentation

Documentation for visualisation interfaces extension module.

## Plotting and Animation

![illustration](assets/jetvis.png)

To visualise the clustered jets as a 3d bar plot (see illustration above)
`Makie.jl` is used. See the `jetsplot` function in `ext/JetVisualisation.jl` and
its documentation for more. There are two worked examples in the `examples`
directory of this package.

The plotting code is a package extension and will load if the one of the `Makie`
modules is loaded in the environment.

The [`animatereco`](@ref) function will animate the reconstruction sequence,
given a `ClusterSequence` object. See the function documentation below for the
many options that can be customised.

## Function Index

```@index
Pages = ["visualisation.md"]
```

## Jet Visualisation Public Interfaces

```@autodocs
Modules = [JetVisualisation]
Order = [:function]
```
