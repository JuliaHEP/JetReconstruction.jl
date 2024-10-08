# EDM4hep Inputs

Extension functionality to read EDM4hep ReconstructedParticles, using the
[EDM4hep.jl](https://github.com/peremato/EDM4hep.jl) package.

## Examples

EDM4hep `ReconstructedParticles` can be used as direct input into jet
reconstruction.

A number of working examples are maintained in the [EDM4hep examples
directory](https://github.com/JuliaHEP/JetReconstruction.jl/tree/main/examples/EDM4hep)
of the package's `examples`.

Here is a snippet that shows the main steps:

```julia
using EDM4hep
using EDM4hep.RootIO
using JetReconstruction

# Change this to something that works on your system
input_file = joinpath("directory", "EDM4hep.root")
reader = RootIO.Reader(input_file)
events = RootIO.get(reader, "events")

evt = events[1]

recps = RootIO.get(reader, evt, "ReconstructedParticles")

cs = jet_reconstruct(recps; algorithm = JetAlgorithm.Durham)
```

## Function Index

```@index
Pages = ["EDM4hep.md"]
```

## EDM4hep Interfaces

```@autodocs
Modules = [EDM4hepJets]
Order = [:function]
```
