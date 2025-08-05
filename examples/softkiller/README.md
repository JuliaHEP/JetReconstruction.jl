# JetReconstruction.jl Examples - SoftKiller

Here are simple examples of SoftKiller that read from two HepMC3 files and
then run the SoftKiller algorithm.

## `softkiller_plots.jl`

This example can be run as:

```sh
julia --project softkiller_plots.jl --pileup-maxevents=100 --eventno=4 --grid-size=0.4 --algorithm=Kt ../../test/data/sk_example_HS.hepmc.zst ../../test/data/sk_example_PU.hepmc.zst
```

This is a simple example of SoftKiller that reads from two HepMC3 files
and displays plots of clustering with SoftKiller and without SoftKiller.

## `softkiller_runtime.jl`

This example can be run as:

```sh
julia --project softkiller_runtime.jl --pileup-maxevents=100 --eventno=4 --grid-size=0.4 --algorithm=Kt ../../test/data/sk_example_HS.hepmc.zst ../../test/data/sk_example_PU.hepmc.zst
```

For both scripts the number of pileup events to add to the hard scatter event
can be controlled with `--pileup-maxevents` and `--pileup-skip`; and the hard
scatter event to work with can be changed with `--eventno`.
