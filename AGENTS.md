# AGENTS.md - Guide for Agentic Coding Assistants

This document provides guidance for agentic coding assistants when working with the **JetReconstruction.jl** package. It outlines the project structure, conventions, and best practices to help AI generate appropriate and consistent code.

Authoritative sources, in order of precedence, if this file disagrees with them:

1. The code in `src/`.
2. `docs/src/contributing.md` - the human contributor guide. **Read this before your first commit**; it carries the PR checklist and the project's development tips.
3. `docs/src/` - the user documentation.

This file exists to get you productive quickly and to warn you about things the
above do not make obvious. It is not a substitute for reading them.

## Package Overview

**JetReconstruction.jl** is a Julia package that implements sequential jet reconstruction (clustering) algorithms for high-energy physics. It reimplements the C++ FastJet package algorithms natively in Julia. It provides additional jet utilities, such as soft killer pileup suppression, substructure operations and Lund plane calculations.

The `README.md` of the project contains basic information.

The `docs/src` directory contains documentation that will be processed by `Documenter.jl`. The build version of the documentation is found [here](https://juliahep.github.io/JetReconstruction.jl/).

The canonical upstream repository is
<https://github.com/JuliaHEP/JetReconstruction.jl>; pull requests are made
against its `main` branch. However, do not commit a branch directly into this
repository, but rather into a fork (usually from the `origin` remote, if the
developer is following the correct forking workflow).

### Nomenclature

Cluster and pseudojet are used interchangeably.

### Sequential Reconstruction

In sequential jet reconstruction final state particles, e.g., measured or simulated in a particle physics detector, are recombined into parent particles that they (putatively) originated from.

The basic steps are outlined here:

- a metric distance is calculated between each pair of clusters in the initial cluster distribution
- for the lowest value of the metric distance, the pair of clusters are merged into a new pseudojet
  - for some algorithms it is possible to do a so-called beam merge, which instead finalises this pseudojet (no more merges are possible), removing it from the list of active clusters
- given the disappearance of two pseudojets and the creation of a new one (or the beam-merge, with one disappearance) metric distances are updated
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

The fixed powers live in `algorithm2power` and the algorithms needing an
explicit `p` in `varpower_algorithms`, both in
`src/AlgorithmStrategyEnums.jl`. Use `is_pp()` / `is_ee()` rather than testing
algorithm identity by hand.

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
  - `src/C_JetReconstruction/` - the C bindings submodule (see *C bindings* below)
- `docs` - package documentation to be processed with `Documenter.jl`
- `test` - package unit and integration tests
  - `test/data` - sample input data and test reference files
- `ext` - additional extension functionality loaded by `Pkg`

Additions:

- `benchmark` - mini-benchmarks for the `JetReconstruction.jl` package (note:
  singular)
- `examples` - examples of using the package, with subdirectories for some
  different functionalities
- `compile` - static compilation of the C library (see *C bindings* below)

### Nested environments

`examples/`, **each of its subdirectories**, `benchmark/`, `docs/` and
`compile/` all carry their own `Project.toml`, each with a
`[sources] JetReconstruction = {path = ".."}` entry (or `"../.."` one level
deeper) pointing back at the package. They are not part of the main project's dependencies. Activate the
right one, and **instantiate it the first time**, or the script will not
resolve:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=examples examples/jetreco.jl --help
```

Each subdirectory is a separate environment needing its own instantiate:

```bash
julia --project=examples/substructure -e 'using Pkg; Pkg.instantiate()'
julia --project=examples/substructure examples/substructure/jet-grooming.jl
```

Skipping the instantiate gives `ArgumentError: Package ArgParse ... is
required but does not seem to be installed`, which looks like a broken script
but is not.

If you add an example needing extra packages, put it in its own
`examples/FEATURE/` subdirectory with its own `Project.toml`.

### Package extensions

`ext/` holds weak-dependency extensions. The exported names are visible without
the trigger package, but do nothing until it is loaded:

| Extension | Trigger package | Provides |
|-----------|-----------------|----------|
| `EDM4hepJets` | `EDM4hep` | EDM4hep input types |
| `JetVisualisation` | `Makie` | `jetsplot`, `animatereco` |
| `JetBenchmarkPlots` | `UnicodePlots` | `plot_trial_times` |

### C bindings

`src/C_JetReconstruction/` exposes a C API, and `compile/` builds it as a
static library (CMake glue in `compile/cmake`, headers in `compile/include`).
`test/test-c-interface.jl` covers it.

**If you change a public signature or the layout of a jet type, check whether
`src/C_JetReconstruction/C_JetReconstruction.jl` and `compile/include` need
updating too.** This is easy to miss and the failure shows up late.

## Development Conventions

### Code Style

- **Formatter**: [JuliaFormatter.jl **v1**](https://github.com/domluna/JuliaFormatter.jl/tree/v1.0.0). v2 formats differently and CI (`.github/workflows/Format.yml`) pins major version 1, so do not let a v2 formatter loose on the tree.
  - Settings come from `.JuliaFormatter.toml`: `style = "sciml"`, `yas_style_nesting = true`, and `*-nb.jl` / `*-pluto.jl` notebooks are excluded. Do not reformat those notebook files.
- **Naming**: Follow [Julia guidelines](https://docs.julialang.org/en/v1/manual/style-guide/): snake_case for functions, PascalCase for types.
- **Type annotations**:
  - Be as concrete as possible with *data members* (use parametrisation if it is needed), do not use abstract data members as these have poor performance
  - Be as general as possible with methods, allowing Julia's natural specialisation to take care of argument types; use abstract struct types to fulfil contracts in the API
- **Docstrings**: Required for all public APIs; add them for non-trivial internal functions too.
- Two further house rules from `docs/src/contributing.md` that are easy to violate:
  - Don't introduce an abstract super type if there will only ever be one concrete sub-type.
  - Don't add accessor methods for purely internal struct fields - access the field directly.

### Julia version floor

`[compat]` sets `julia = "1.10"`, and CI runs `lts`, `1` and `nightly` on both
x64 and aarch64. **Do not use language features newer than 1.10.** A
`Downgrade compat` workflow additionally tests against the lowest permitted
version of every dependency, so do not rely on behaviour newer than the
declared `[compat]` floor of a dependency either.

### Adding a dependency

`test/test-aqua.jl` runs `Aqua.test_all`. In practice this means:

- **A new entry in `[deps]` must have a matching `[compat]` entry**, or CI
  fails. Same for `[extras]` used by tests.
- No method ambiguities, no type piracy, no stale dependencies.

### Type Parameters

The package makes extensive use of type parameters for performance. The one you
will meet constantly is:

- `T <: FourMomentum` - the *internal jet type*, concretely `PseudoJet` (pp) or
  `EEJet` (e+e-). `FourMomentum` is declared in `src/CommonJet.jl`. Note this is
  the type inputs are **converted into**, not the user's input type;
  `ClusterSequence{T}` and most recombination methods are parametrised on it.

User input particles are a separate matter - see *Input Particle Types* below.
Other type parameters are local to their own file (for example `Surrounding{N}`
in `src/TiledAlgoLLStructs.jl`, where `N` is a count of neighbouring tiles) and
carry no package-wide meaning.

### Performance Considerations

- Use `@inbounds` and `@simd` where appropriate - this is the dominant idiom in the hot loops of `src/PlainAlgo.jl`, `src/EEAlgorithm.jl` and `src/TiledAlgoLL.jl`.
- Prefer `StructArrays.jl` for arrays-of-structs to keep hot loops columnar
- Avoid dynamic dispatch in inner loops; keep functions type-stable
- `LoopVectorization.jl` is a dependency but is used in exactly one place where it makes a significant difference (`@turbo` in `src/Utils.jl:fast_findmin()`). Do not use it by default - match the surrounding code and measure performance first.
- Profile before optimising. See *Benchmarking and profiling* below.

### Tests

All new features should be supported by tests. The pattern used is that feature
`foo` should be tested in the specific test file `test/test-foo.jl`.

**You must then add `include("test-foo.jl")` to `test/runtests.jl` yourself** -
it is a hand-maintained list, nothing globs the directory.

Note the use of `test/common.jl`, a guard wrapper that lets sub-tests run
independently; it includes the real fixtures in `test/_common.jl`. Every test
file should start with `include("common.jl")`.

#### Reference Data

Any reference data should live in `test/data`. Zstd compression is preferred;
`JetReconstruction.open_with_stream()` (in `src/Utils.jl`) makes the
decompression transparent.

#### Running tests

The full suite:

```bash
julia --project test/runtests.jl
```

or from the REPL, `] test JetReconstruction`.

A single test file standalone:

```bash
julia --project test/test-foo.jl
```

**This requires `TestEnv.jl` to be installed in your *default* environment** -
`test/common.jl` calls `TestEnv.activate()` when run as a script. Without it you
get an `ArgumentError` about the `TestEnv` package, which has nothing to do with
the test itself. Install once with:

```bash
julia -e 'using Pkg; Pkg.activate(); Pkg.add("TestEnv")'
```

### Documentation

The API reference pages (`docs/src/lib/public.md` and `internal.md`) are
`@autodocs` blocks, so **a newly exported symbol needs no manual edit there** -
a good docstring is enough. A new *narrative* page does have to be added to the
`pages` list in `docs/make.jl`.

To build:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=docs docs/make.jl
```

This pulls in CairoMakie, EDM4hep and UnicodePlots and is slow; do not build the
docs just to check a docstring renders.

### Benchmarking and profiling

`BenchmarkTools` is **not** a dependency of the main project - `@btime` will not
resolve under `--project`. It lives in `benchmark/Project.toml`:

```bash
julia --project=benchmark benchmark/runbenchmarks.jl --verbose benchmark_results.json
```

For profiling, `examples/instrumented-jetreco.jl` is the worked example, and
`examples/Project.toml` carries `Profile`, `FlameGraphs` and
`StatProfilerHTML`.

### Git and pull requests

- Develop against `main` (see `docs/src/contributing.md` for the exception on
  `release-X` branches for bug fixes).
- If asked to push to remote, this should be the developer's `origin` fork of the upstream repository
- Open an issue to discuss anything large before implementing it.
- Run the formatter and the full test suite before proposing a PR.
- Rebase on `main` if review takes a while and `main` had advanced.

## Common Patterns

These are quick orientation snippets. `docs/src/` is the maintained reference;
if it disagrees with what follows, believe the docs and please fix this file.

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
- `EDM4hep.ReconstructedParticle` (requires the `EDM4hep` extension)

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

The defaults are `preprocess_escheme` / `addjets_escheme` in `src/CommonJet.jl`.

### Using Named Recombination Schemes

```julia
using JetReconstruction: RecombinationScheme, RecombinationMethods

myscheme = RecombinationMethods[RecombinationScheme.PtScheme]
cs = jet_reconstruct(particles; R=1.0, algorithm=JetAlgorithm.AntiKt, myscheme...)
```

Available schemes: `EScheme`, `ESchemeRaw`, `PtScheme`, `Pt2Scheme`

### Jet Substructure

**All substructure entry points are `PseudoJet`-only**: they are typed
`(jet::PseudoJet, clusterseq::ClusterSequence{PseudoJet})`. They will not accept
an `EEJet` cluster sequence.

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

### SoftKiller and Lund plane

Not summarised here - see `docs/src/softkiller.md` and `docs/src/lundplane.md`,
with worked examples under `examples/softkiller/` and `examples/lundplane/`.

## Common Pitfalls

1. **Aqua**: a new dependency without a `[compat]` entry breaks CI
2. **Nested environments**: example, benchmark, docs and compile scripts need their own `--project=...`
3. **Standalone tests**: need `TestEnv.jl` in your default environment
4. **New test files**: must be added to `test/runtests.jl` by hand
5. **Type stability**: ensure functions are type-stable, especially in hot paths
6. **Memory allocation**: minimize allocations in inner loops
7. **Algorithm parameters**: remember that Durham ignores R, Valencia requires β and γ
8. **Strategy selection**: N2Plain is for low particle counts, N2Tiled for high counts
9. **Particle types**: pp uses `PseudoJet`, e+e- uses `EEJet`; substructure is `PseudoJet`-only
10. **Coordinate systems**: be consistent with rapidity vs pseudorapidity
11. **C bindings**: a public API change may require a matching change in `src/C_JetReconstruction/` and `compile/include`
12. **Julia 1.10 floor**: no newer language features

## Useful Commands

```bash
# Run tests
julia --project test/runtests.jl

# Run a specific test standalone (needs TestEnv.jl in the default environment)
julia --project test/test-feature.jl

# Format code (JuliaFormatter v1 - bootstrap the shared env once, see below)
julia --project=@juliaformatter -e 'using JuliaFormatter; format(".")'

# One-time bootstrap of that formatter environment
julia --project=@juliaformatter -e 'using Pkg; Pkg.add("JuliaFormatter"); Pkg.compat("JuliaFormatter", "1.0"); Pkg.update()'

# Build documentation (slow)
julia --project=docs docs/make.jl

# Run an example
julia --project=examples examples/jetreco.jl --help

# Run the benchmarks
julia --project=benchmark benchmark/runbenchmarks.jl --verbose benchmark_results.json
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
| EDM4hep integration | `examples/EDM4hep/` (SimpleRecoEDM4hep.jl, EDM4hepJets.jl, QuickLookEDM4hep.jl) |
| Jet substructure (grooming, tagging) | `examples/substructure/` (jet-grooming.jl, jet-tagging.jl) |
| Constituent access | `examples/constituents/` (jetreco-constituents.jl, jetreco-constituents-nb.jl) |
| Lund plane | `examples/lundplane/` (lund-jet-generation.jl, lund-plane-visualisation.jl) |
| SoftKiller pileup mitigation | `examples/softkiller/` (softkiller_plots.jl, softkiller_runtime.jl) |
| C API usage | `test/test-c-interface.jl`, `compile/downstream` |

## References

- FastJet: <https://fastjet.fr/>
- ArXiv papers: hep-ph/0512210, arXiv:1111.6097, 1404.4294
- "Fast Jet Finding in Julia", EPJ Web Conf. 337 (2025) 01067, <https://arxiv.org/abs/2503.08146>
- JuliaHEP ecosystem: <https://github.com/JuliaHEP>
- Documentation: <https://juliahep.github.io/JetReconstruction.jl/>

## When in Doubt

0. Ask the human, with an explanation of the issue and the options for moving forward
1. Check existing implementations in `src/` for patterns
2. Look at test files for usage examples
3. Consult `docs/src/contributing.md` and the documentation at <https://juliahep.github.io/JetReconstruction.jl/>
4. Follow Julia conventions and the Julia Style Guide
5. Profile before optimizing
