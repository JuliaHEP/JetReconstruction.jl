# AGENTS.md - Guide for Agentic Coding Assistants

This document provides guidance for agentic coding assistants when working with the **JetReconstruction.jl** package. It outlines the project structure, conventions, and best practices to help AI generate appropriate and consistent code.

## Package Overview

**JetReconstruction.jl** is a Julia package that implements sequential jet reconstruction (clustering) algorithms for high-energy physics. It reimplements the C++ FastJet package algorithms natively in Julia. It provides additional jet utilities, such as soft killer pileup suppression, substructure operations and Lund plane calculations.

The @README.md of the project contains basic information.

The `docs/src` directory contains documentation that will be processed by `Documenter.jl`. The build version of the documentation is found [here](https://juliahep.github.io/JetReconstruction.jl/).

### Nomenclature

Cluster and pseudojet are used interchangeably.

### Sequential Reconstruction

In sequential jet reconstruction final state particles, e.g., measured or simulated in a particle physics detector, are recombined into parent particles that they (putatively) originated from.

The basic steps are outlined here:

- a metric distance is calculated between each pair of clusters in the initial cluster distribution
- for the lowest value of the metric distance, the pair of clusters are merged into a new pseudojet
  - for some algorithms it is possible to do a so-called beam merge, which instead finalises this pseudojet (no more merges are possible), removing it from the list of active clusters
- given the disappearance of two pesudojets and the creation of a new one (or the beam-merge, with one disappearance) metric distances are updated
- the loop then continues, with the next lowest distance pseudojets being merged, until no more merges are possible and all initial clusters have been processed

Note that special attention is paid in the implementation to minimise the calculations performed, so as to increase speed (see `Reconstruction Strategies` below.)

### Key Algorithms

Algorithms express different mathematical formula used to calculate the metric distance. e.g., they may set a different power value (p) on the momentum scaling used when calculating the distance.

The other major difference is that in pp reconstruction the (rapidity, phi) space is used; in the e+e- reconstruction the (theta, phi) space is used for the geometric distance component.

| Algorithm | Julia EnumX Type | Notes |
|-----------|------|-------|
| Anti-kₜ | `JetAlgorithm.AntiKt` | Default for pp collisions, p=-1, R value controls radius of maximum merger distance |
| Cambridge/Aachen | `JetAlgorithm.CA` | pp collisions, p=0, R used |
| kₜ (inclusive) | `JetAlgorithm.Kt` | pp collisions, p=1, R used |
| Generalised kₜ (pp) | `JetAlgorithm.GenKt` | pp collisions, arbitrary `p` parameter, R used |
| Durham | `JetAlgorithm.Durham` | e+e- collisions, p=1, no R value |
| Generalised kₜ (e+e-) | `JetAlgorithm.EEKt` | e+e- collisions, arbitrary `p` parameter, R used |
| Valencia | `JetAlgorithm.Valencia` | e+e- collisions, requires `p` (β) and `γ` |

### Reconstruction Strategies

For *pp algorithms only* there are strategy options that can be selected. Usually `RecoStrategy.Best` is used.

| Strategy | Type | Use Case |
|----------|------|----------|
| Best | `RecoStrategy.Best` | Auto-selects based on particle density |
| N2Plain | `RecoStrategy.N2Plain` | Global matching, good for low N |
| N2Tiled | `RecoStrategy.N2Tiled` | Tiled search, good for high N |

The e+e- algorithms and N2Plain share the same basic structure, in that all particles are measured for nearest neighbour distances. In the N2Tiled case there is bookkeeping of particles into tiles, based on R, that reduce the number of possible matches of particles - this means the algorithm scales much better for higher cluster densities.

## Project Structure

This project follows standard Julia conventions on layout, with a few additions:

- `src` - package source files, with the main package entry point `src/JetReconstruction.jl`
- `docs` - package documentation to be processes with `Documenter.jl`
- `test` - package unit and integration tests
  - `test/data` - sample input data and test reference files
- `ext` - additional extension functionality loaded by `Pkg`

Additions:

- `benchmarks` - mini-benchmarks for the `JetReconstruction.jl` package
- `examples` - examples of using the package, with subdirectories for some different functionalities (note, has it's own `Package.toml` for packages that are not required in the main package itself)

## Development Conventions

### Code Style

- **Formatter**: Use [JuliaFormatter.jl v1](https://github.com/domluna/JuliaFormatter.jl/tree/v1.0.0)
- **Naming**: Follow [Julia guidelines](https://docs.julialang.org/en/v1/manual/documentation/#Writing-Documentation): snake_case for functions, PascalCase for types.
- **Type annotations**: 
  - Be as concrete as possible with *data members* (use parametrisation if it is needed), do not use abstract data members as these have poor performance
  - Be as general as possible with methods, allowing Julia's natural specialisation to take care of argument types; use abstract struct types to fulfil contracts in the API
- **Docstrings**: Required for all public APIs

### Type Parameters

The package makes extensive use of type parameters for performance. Key type variables:

- `T` - Input particle type (must implement LorentzVectorBase interface)
- `N` - Number of dimensions (typically 4 for spacetime)
- `S` - Storage type for LorentzVectorHEP

### Performance Considerations

- Use `@inbounds` and `@simd` where appropriate
- Prefer `StructArrays.jl` for arrays-of-structs to keep hot loops columnar
- Use `LoopVectorization.jl` for hot loops
- Avoid dynamic dispatch in inner loops
- Profile with `@time`, `@btime`, and `--project` profiling tools

### Tests

All new features should be supported by tests. The pattern used is that feature
`foo` should be tested in the specific test file `test/test-foo.jl`, which will
be included in the main test suite entry point `test/runtests.jl`.

Note the use of `test/common.jl` which allows sub-tests to run independently.

#### Reference Data

Any reference data should live in `test/data`.

#### Running tests

For a specific feature:

```julia
julia --project test/test-foo.jl
```

or for the full suite:

```
julia
julia --project test/runtests.jl
```

## Common Patterns

### Basic Jet Reconstruction

```julia
using JetReconstruction

# For pp collisions
particles = [...]  # Vector of PseudoJet or any LorentzVectorBase type
cs = jet_reconstruct(particles; algorithm=JetAlgorithm.AntiKt, R=1.0)
jets = inclusive_jets(cs; ptmin=5.0)

# For e+e- collisions
cs = jet_reconstruct(particles; algorithm=JetAlgorithm.Durham)
jets = exclusive_jets(cs; njets=4)
```

### Input Particle Types

Any type implementing the `LorentzVectorBase.jl` interface will work, i.e., types `T` that provide `LorentzVectorBase.coordinate_system(::T)`.

Supported types include:

- `PseudoJet` (for pp)
- `EEJet` (for e+e-)
- `LorentzVectorHEP.LorentzVector`
- `EDM4hep.ReconstructedParticle`

### Custom Recombination Schemes

```julia
# Define preprocessing function
function my_preprocess(jet::T, ::Type{OutputT}; cluster_hist_index) where {T, OutputT}
    OutputT(px(jet), py(jet), pz(jet), energy(jet); cluster_hist_index=cluster_hist_index)
end

# Define recombination function  
function my_recombine(jet1::T, jet2::T; cluster_hist_index::Int) where {T}
    T(px(jet1) + px(jet2),
      py(jet1) + py(jet2),
      pz(jet1) + pz(jet2),
      energy(jet1) + energy(jet2);
      cluster_hist_index=cluster_hist_index)
end

# Use in reconstruction
cs = jet_reconstruct(particles; 
                      algorithm=JetAlgorithm.AntiKt, 
                      R=1.0,
                      preprocess=my_preprocess,
                      recombine=my_recombine)
```

### Using Named Recombination Schemes

```julia
using JetReconstruction: RecombinationScheme, RecombinationMethods

myscheme = RecombinationMethods[RecombinationScheme.PtScheme]
cs = jet_reconstruct(particles; R=1.0, algorithm=JetAlgorithm.AntiKt, myscheme...)
```

Available schemes: `EScheme`, `ESchemeRaw`, `PtScheme`, `Pt2Scheme`

### Jet Substructure

```julia
# Mass Drop
mass_drop(jet, clusterseq; mu=0.67, y=0.09)

# Soft Drop
soft_drop(jet, clusterseq; zcut=0.1, beta=2.0, radius=1.0)

# Filtering
jet_filtering(jet, clusterseq; radius=0.3, hardest_jets=3)

# Trimming
jet_trimming(jet, clusterseq; radius=0.3, fraction=0.3, recluster_method=JetAlgorithm.CA)
```

### Accessing Cluster Sequence Information

```julia
# Get inclusive jets with pT cut
jets = inclusive_jets(cs; ptmin=10.0)

# Get exclusive jets
jets = exclusive_jets(cs; dcut=1.0)  # or njets=3

# Get jet constituents (indexes into original particles)
indexes = constituent_indexes(jet, cs)

# Get jet constituents (actual particle objects)
constituents_list = constituents(jet, cs)

# Get parent jets
parents = parent_jets(jet, cs)  # Returns Tuple or nothing
```

### Documentation Generation

The docs use Documenter.jl. To build:

```bash
julia --project=docs docs/make.jl
```

## Common Pitfalls

1. **Type stability**: Ensure functions are type-stable, especially in hot paths
2. **Memory allocation**: Minimize allocations in inner loops
3. **Algorithm parameters**: Remember that Durham ignores R, Valencia requires β and γ
4. **Strategy selection**: N2Plain is for low particle counts, N2Tiled for high counts
5. **Particle types**: pp uses `PseudoJet`, e+e- uses `EEJet`
6. **Coordinate systems**: Be consistent with rapidity vs pseudorapidity

## Useful Commands

```bash
# Format code
julia --project=@juliaformatter -e 'using JuliaFormatter; format(".")'

# Run tests
julia --project test/runtests.jl

# Run specific test
julia --project test/test-feature.jl

# Build documentation
julia --project=docs docs/make.jl

# Profile code
julia --project -e 'using Profile; @profile include("script.jl")'

# Memory profiling
julia --project -e '@time include("script.jl")'
```

## Where to Look for Examples

| Task | Location |
|------|----------|
| Basic reconstruction | `examples/jetreco.jl` |
| Performance profiling | `examples/instrumented-jetreco.jl` |
| Command-line options parsing | `examples/parse-options.jl` |
| 3D jet visualisation | `examples/visualisation/visualise-jets.jl` |
| Reconstruction animation | `examples/visualisation/animate-reconstruction.jl` |
| Jupyter notebook (visualisation) | `examples/visualisation/visualise-jets.ipynb` |
| Pluto notebook (visualisation) | `examples/visualisation/visualise-jets-nb.jl` |
| EDM4hep integration | `examples/EDM4hep/` (SimpleRecoEDM4hep.jl, EDM4hepJets.jl) |
| Jet substructure (grooming, tagging) | `examples/substructure/` (jet-grooming.jl, jet-tagging.jl) |
| Constituent access | `examples/constituents/` (jetreco-constituents.jl, jetreco-constituents-nb.jl) |
| Lund plane | `examples/lundplane/` (lund-jet-generation.jl, lund-plane-visualisation.jl) |
| SoftKiller pileup mitigation | `examples/softkiller/` (softkiller_plots.jl, softkiller_runtime.jl) |

## References

- FastJet: <https://fastjet.fr/>
- ArXiv papers: hep-ph/0512210, arXiv:1111.6097, 1404.4294
- JuliaHEP ecosystem: <https://github.com/JuliaHEP>
- Documentation: <https://juliahep.github.io/JetReconstruction.jl/>

## When in Doubt

0. Ask the human, with an explanation of the issue and the options for moving forward
1. Check existing implementations in `src/` for patterns
2. Look at test files for usage examples
3. Consult the documentation at <https://juliahep.github.io/JetReconstruction.jl/>
4. Follow Julia conventions and the Julia Style Guide
5. Profile before optimizing
