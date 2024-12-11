# Jet Reconstruction Examples

The Jet Reconstruction package has a number of example files that show how to
usage. These are in the `examples` subdirectory of the package and can be
browsed directly on
[GitHub](https://github.com/JuliaHEP/JetReconstruction.jl/tree/main/examples).

*Note:* because of extra dependencies in these scripts, one must use the
`Project.toml` file in the `examples` directory.

## Standalone Examples

### `jetreco.jl`

This is a basic jet reconstruction example that shows how to call the package to
perform a jet reconstruction, with different algorithms and (optionally)
strategy, producing exclusive and inclusive jet selections.

```sh
julia --project=examples examples/jetreco.jl --algorithm=AntiKt test/data/events.pp13TeV.hepmc3.gz
...
julia --project=examples examples/jetreco.jl --algorithm=Durham test/data/events.eeH.hepmc3.gz
...
julia --project=examples examples/jetreco.jl --maxevents=10 --strategy=N2Plain --algorithm=Kt --exclusive-njets=3 test/data/events.pp13TeV.hepmc3.gz
...
```

There are options to explicitly set the algorithm (use `--help` to see these).

### `instrumented-jetreco.jl`

This is a more sophisticated example that allows performance measurements to be
made of the reconstruction, as well as profiling (flamegraphs and memory
profiling). Use the `--help` option to see usage. e.g., to extract timing
performance for the AntiKt algorithm using the tiled strategy:

```sh
julia --project instrumented-jetreco.jl -S N2Tiled -A AntiKt --nsamples 100 ../test/data/events.hepmc3
```

### `visualise-jets.jl`

This script will produce a PNG/PDF showing the results of a jet reconstruction.
This is a 3D plot where all the initial energy deposits are visualised, with
colours that indicate in which final cluster the deposit ended up in.

### `visualise-jets.ipynb`

Similar to `visualise-jets.jl` this notebook will produce a visualisation of jet
reconstruction in the browser. This is a 3D plot where all the initial energy
deposits are visualised, with colours that indicate in which final cluster the
deposit ended up in.

### `animate-reconstruction.jl`

Performs jet reconstruction and then produces and animation of the process,
showing how the jets merge from their different constituents.

## EDM4hep

The `examples/EDM4hep` folder contains examples of using EDM4hep reconstructed
particles as input to jet reconstruction. See the specific `README.md` file in
that directory as well as [EDM4hep Inputs](@ref).

## Jet Constituents

The `examples/constituents` folder shows an example of the two mechanisms to
retrieve jet constituents.
