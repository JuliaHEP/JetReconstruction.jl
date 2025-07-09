# JetReconstruction.jl Examples - SoftKiller

## `softkiller_plots.jl`

This example can be run as:

```sh
julia --project softkiller_plots.jl  --maxevents=100 --grid-size=0.4  --algorithm=Kt --pileup-file=../../test/data/sk_example_PU.hepmc.zst --hard-file=../../test/data/sk_example_HS.hepmc.zst
```

This is a simple example of SoftKiller that reads from two HepMC3 files
and displays plots of clustering without SoftKiller and with SoftKiller.

## `softkiller_runtime.jl`

This example can be run as:

```sh
julia --project softkiller_runtime.jl  --maxevents=100 --grid-size=0.4  --algorithm=Kt --pileup-file=../../test/data/sk_example_PU.hepmc.zst --hard-file=../../test/data/sk_example_HS.hepmc.zst
```

This is a simple example of SoftKiller that reads from two HepMC3 files
and displays the runtime of the SoftKiller algorithm.
