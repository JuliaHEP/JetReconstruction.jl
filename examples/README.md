# JetReconstruction.jl Examples

This directory and subdirectories have a number of example files that show how
to used the `JetReconstruction.jl` package.

Because of extra dependencies in these scripts, one must use the `Project.toml`
file in this directory.

Some features are demonstrated in their own subdirectories, in which case use
the `Project.toml` file in that folder.

## `jetreco.jl`

This is a basic jet reconstruction example that shows how to call the package to
perform a jet reconstruction, with different algorithms and (optionally)
strategy, producing exclusive and inclusive jet selections.

## `instrumented-jetreco.jl`

This is a more sophisticated example that allows performance measurements to be
made of the reconstruction, as well as profiling (flamegraphs and memory
profiling).
