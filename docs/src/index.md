# Jet Reconstruction

This package implements sequential Jet Reconstruction (clustering) algorithms,
which are used in high-energy physics as part of event reconstruction for ``pp``
and ``e^+e^-`` colliders.

## Algorithms

Algorithms used are based on the C++ FastJet package (<https://fastjet.fr>,
[hep-ph/0512210](https://arxiv.org/abs/hep-ph/0512210),
[arXiv:1111.6097](https://arxiv.org/abs/1111.6097)), reimplemented natively in
Julia.

The algorithms include anti-``{k}_\text{T}``, Cambridge/Aachen, inclusive
``k_\text{T}``, generalised ``k_\text{T}`` for ``pp`` events; and the Durham
algorithm and generalised ``k_\text{T}`` for ``e^+e^-``. The Valencia algorithm
for ``e^+e^-``, based on [1404.4294](https://arxiv.org/abs/1404.4294) is also
supported.

## Reconstruction Interface

The main interface for reconstruction is [`jet_reconstruct`](@ref), called as, e.g.,

```julia
jet_reconstruct(particles; algorithm = JetAlgorithm.AntiKt, R = 1.0)
```

or with some of the optional arguments,

```julia
jet_reconstruct(particles; algorithm = JetAlgorithm.GenKt, R = 0.4, 
                p = 0.5, recombine = addjets, strategy = RecoStrategy.Best)
```

Where `particles` is a collection of 4-vector objects (see [Input Particle
Types](@ref)) to reconstruct and the algorithm is given explicitly.

For the case of generalised ``k_T`` (for ``pp`` and ``e^+e^-``) both the
algorithm (`GenKt`, `EEKt`) and `p` are needed.

For the Valencia algorithm, as well as `R`, the β (equivevent to the existing
power p of the algorithm) and γ (angular exponent parameter used in the beam
distance) are needed.

The `R` value determines the cone size; in the case of the Durham algorithm the
`R` value is ignored.

For a discussion of the `recombine` function, see [Jet Recombination](@ref).

The object returned is a [`ClusterSequence`](@ref), which internally tracks all
merge steps and is used for [Inclusive and Exclusive Selections](@ref).

### Algorithm Types

Each known algorithm is referenced using a `JetAlgorithm` scoped enum value.

| Algorithm | Type name | Notes |
|-----------|-----------|-------|
| anti-``{k}_\text{T}`` | `JetAlgorithm.AntiKt` | Implies `p=-1` |
| Cambridge/Aachen | `JetAlgorithm.CA` | Implies `p=0` |
| inclusive ``k_\text{T}`` | `JetAlgorithm.Kt` | Implies `p=1` |
| generalised ``k_\text{T}`` | `JetAlgorithm.GenKt` | For $pp$, value of `p` must also be specified |
| ``e^+e-`` ``k_\text{T}`` / Durham | `JetAlgorithm.Durham` | `R` value ignored and can be omitted |
| generalised ``e^+e-`` ``k_\text{T}`` | `JetAlgorithm.EEKt` | For ``e^+e^-``, value of `p` must also be specified |
| Valencia | `JetAlgorithm.Valencia` | For ``e^+e^-``, values of `p` (β) and `γ` must be specified |

### Strategy

Generally one does not need to manually specify a strategy, but [Algorithm
Strategy](@ref) describes how to do this, if desired.

## Inclusive and Exclusive Selections

To obtain final jets both inclusive (``p_T`` cut) and exclusive (``n_{jets}`` or
``d_{ij}`` cut) selections are supported:

- [`inclusive_jets(clusterseq::ClusterSequence, ptmin = 0.0)`](@ref)
- [`exclusive_jets(clusterseq::ClusterSequence; dcut = nothing, njets = nothing)`](@ref)

(For `exclusive_jets` either `dcut` or `njets` is needed, but not both.)

### Sorting

Sorting vectors is trivial in Julia, no special sorting methods are provided. As
an example, to sort exclusive jets of ``>5.0`` (usually GeV, depending on your
EDM) from highest energy to lowest:

```julia
sorted_jets = sort!(inclusive_jets(cs::ClusterSequence; ptmin=5.0), 
  by=JetReconstruction.energy, rev=true)
```

## Jet Constituents and Jet Parents

There are two ways to retrieve jet constituents. The first way is just to
retrieve the *indexes* of the constituent jets. These indexes refer to the
original collection of particles passed in to the reconstruction.

- [`constituent_indexes`](@ref)

The alternative it to retrieve the actual jets from the reconstruction sequence.
In this case the returned array contains references to the jet objects (of type
`T`) used internally in the reconstruction.

- [`constituents`](@ref)

Note that in both these cases the cluster sequence object from the
reconstruction is required (to avoid circular dependencies and improve memory
management reconstructed jets do not contain a link back to their cluster
sequence).

To retrieve a jet's parents:

- [`parent_jets`](@ref)

This will return a tuple of the target jet's parents, or `nothing` when one or
both parents are missing (the only case when a jet has one parent is when it
undergoes a *beam merge* step).

## References


The current recommended reference for JetReconstruction.jl is:

```bibtex
@article{refId0,
  author = {{Stewart, Graeme Andrew} and {Ganguly, Sanmay} and {Ghosh, Sattwamo} and {Gras, Philippe} and {Krasnopolski, Atell}},
  title = {Fast Jet Finding in Julia},
  DOI= "10.1051/epjconf/202533701067",
  url= "https://doi.org/10.1051/epjconf/202533701067",
  journal = {EPJ Web Conf.},
  year = 2025,
  volume = 337,
  pages = "01067",
}
```

Also available as [arXiv:2503.08146](https://arxiv.org/abs/2503.08146).

### Other Articles

- CHEP2023, *Polyglot Jet Finding*: [arXiv:2309.17309](https://arxiv.org/abs/2309.17309), [10.1051/epjconf/202429505017](https://doi.org/10.1051/epjconf)

## Community

We welcome contributions to the project (see [Contributing to
JetReconstruction](@ref)). Please follow our [Code of
Conduct](https://github.com/JuliaHEP/JetReconstruction.jl/blob/main/CODE_OF_CONDUCT.md),
thanks!

## Authors and Copyright

Code in this package is authored by:

- Atell Krasnopolski <delta_atell@protonmail.com>
- Graeme A Stewart <graeme.andrew.stewart@desy.de>
- Philippe Gras <philippe.gras@cern.ch>

and is Copyright 2022-2025 The Authors, CERN.

The code is under the MIT License.
