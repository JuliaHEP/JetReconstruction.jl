# Jet Reconstruction from EDM4hep inputs

The `EDM4hepJets.jl` file shows an example of how to reconstruct jets from and
[EDM4hep](https://edm4hep.web.cern.ch) input file. The reconstruction is set to
use the Durham algorithm and to output 2 exclusive jets.

To use this example you need to give an EDM4hep input file as an argument, e.g.,

```sh
julia --project EDM4hepJets.jl [--maxevents MAX_EVENTS] path/to/input.root
```

The optional argument `--maxevents` controls how many events are read from the
input file (default is just 1).

There is also a very simple `SimpleRecoEDM4hep.jl` that can be used in the REPL
for tests and development. Here one has to update the file path explicitly to a
valid local EDM4hep file.
