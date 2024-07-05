# Jet Reconstruction Examples

The Jet Reconstruction package has a number of example files that show how to
usage.

*Note:* because of extra dependencies in these scripts, one must use the
`Project.toml` file in the `examples` directory.

## `jetreco.jl`

This is a basic jet reconstruction example that shows how to call the package to
perform a jet reconstruction, with different algorithms and (optionally)
strategy, producing exclusive and inclusive jet selections.

```sh
julia --project=. examples/jetreco.jl --maxevents=100 --strategy=N2Plain test/data/events.hepmc3
...
julia --project=. examples/jetreco.jl --maxevents=100 --strategy=N2Tiled test/data/events.hepmc3
...
```

There are options to explicitly set the algorithm (use `--help` to see these).

## `instrumented-jetreco.jl`

This is a more sophisticated example that allows performance measurements to be
made of the reconstruction, as well as profiling (flamegraphs and memory
profiling). Use the `--help` option to see usage. e.g., to extract timing
performance for the AntiKt algorithm using the tiled strategy:

```sh
julia --project instrumented-jetreco.jl -S N2Tiled -A AntiKt --nsamples 100 ../test/data/events.hepmc3
```

## `visualise-jets.jl`

This script will produce a PNG/PDF showing the results of a jet reconstruction.
This is a 3D plot where all the initial energy deposits are visualised, with
colours that indicate in which final cluster the deposit ended up in.

## `visualise-jets.ipynb`

Similar to `visualise-jets.jl` this notebook will produce a visualisation of jet
reconstruction in the browser. This is a 3D plot where all the initial energy
deposits are visualised, with colours that indicate in which final cluster the
deposit ended up in.

## `animate-reconstruction.jl`

Performs jet reconstruction and then produces and animation of the process,
showing how the jets merge from their different constituents.
