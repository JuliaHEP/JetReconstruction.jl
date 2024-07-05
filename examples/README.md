# JetReconstruction.jl Examples

This directory has a number of example files that show how to used the
`JetReconstruction.jl` package.

Because of extra dependencies in these scripts, one must use the `Project.toml`
file in this directory.

## `jetreco.jl`

This is a basic jet reconstruction example that shows how to call the package to
perform a jet reconstruction, with different algorithms and (optionally)
strategy, producing exclusive and inclusive jet selections.

## `instrumented-jetreco.jl`

This is a more sophisticated example that allows performance measurements to be
made of the reconstruction, as well as profiling (flamegraphs and memory
profiling).

## `visualise-jets.jl`

This script will produce a PNG/PDF showing the results of a jet reconstruction.
This is a 3D plot where all the initial energy deposits are visualised, with
colours that indicate in which final cluster the deposit ended up in.

## `visualise-jets.ipynb`

This notebook will produce a visualisation of jet reconstruction in the browser.
This is a 3D plot where all the initial energy deposits are visualised, with
colours that indicate in which final cluster the deposit ended up in.

## `animate-reconstruction.jl`

Performs jet reconstruction and then produces and animation of the process,
showing how the jets merge from their different constituents.
