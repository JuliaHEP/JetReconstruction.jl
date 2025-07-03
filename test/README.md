# Tests for JetReconstruction

This directory contains the unit tests for the JetReconstruction package.

To run all of the tests do:

```sh
julia --project test/runtests.jl
```

or

```sh
julia --project
julia> ]
(JetReconstruction) pkg> test
...
```

Tests are factorised into individual files as `test-*.jl`. By use of `common.jl` all of these tests can be run standalone, provided that `TestEnv.jl` is installed in the default environment, e.g.,

```sh
julia --project test/test-pp-reconstruction.jl
```

## Data Files

In `data`:

- `events.pp13TeV.hepmc3.zst` - 100 Pythia 13 TeV events in HepMC3 format (compressed)
- `events.eeH.hepmc3.zst` - 100 Pythia $e^+e^- \rightarrow H$ events in HepMC3 format (compressed)

- `hepmc32summary.jl` script that will calculate the average particle density of
  a HepMC3 input file and plot an ASCII histogram of the distribution.

Other files are reference outputs used by tests.
