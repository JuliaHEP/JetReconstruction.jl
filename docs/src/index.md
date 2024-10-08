```@raw html
---
layout: home

hero:
  name: "JetReconstruction.jl"
  tagline:  "Jet reconstruction (reclustering) with Julia"
  image:
    src: /logo.png
    alt: DocumenterVitepress
  actions:
    - theme: brand
      text: Examples
      link: /examples
    - theme: alt
      text: Public APIs
      link: /lib/public

---
```

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
algorithm and generalised ``k_\text{T}`` for ``e^+e^-``.

## Reconstruction Interface

The main interface for reconstruction is [`jet_reconstruct`](@ref), called as, e.g.,

```julia
jet_reconstruct(particles; algorithm = JetAlgorithm.AntiKt, R = 1.0)
```

or with some of the optional arguments,

```julia
jet_reconstruct(particles; algorithm = JetAlgorithm.GenKt, R = 0.4, p = 0.5, recombine = +, strategy = RecoStrategy.Best)
```

Where `particles` is a collection of 4-vector objects to reconstruct and the
algorithm is either given explicitly or implied by the power value. For the case
of generalised ``k_T`` (for ``pp`` and ``e^+e^-``) both the algorithm
(`JetAlgorithm.GenKt`) and `p` are needed. For the case of the Durham algorithm
the `R` value is ignored.

The object returned is a [`ClusterSequence`](@ref), which internally tracks all
merge steps.

### Input Particle Types

For the `particles` input to the reconstruction any one dimensional
`AbstractArray{T, 1}` can be used, where the type `T` has to implement methods
to extract the 4-vector components, viz, the following are required:

- `JetReconstuction.px(particle::T)`
- `JetReconstuction.py(particle::T)`
- `JetReconstuction.pz(particle::T)`
- `JetReconstuction.energy(particle::T)`

Currently built-in supported types are
[`LorentzVectorHEP`](https://github.com/JuliaHEP/LorentzVectorHEP.jl), the
`PseudoJet` and `EEjet`s from this package, and `ReconstructedParticles` from
[EDM4hep](https://github.com/peremato/EDM4hep.jl).

If you require support for a different input collection type then ensure you
define the `px()`, etc. methods *for your specific type* and *in the
`JetReconstruction` package*. This use of what might be considered type piracy
is blessed as long as you are en *end user* of the jet reconstruction package.

If your type is used in several places or by different users, please consider
writing a package extension that will support your type, following the model for
EDM4hep in `ext/EDM4hepJets.jl`.

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

### ``pp`` Algorithms

For the three ``pp`` algorithms with fixed `p` values, the `p` value can be
given instead of the algorithm name. However, this should be considered
*deprecated* and will be removed in a future release.

## Strategy

For the ``pp`` algorithms three strategies are available for the different
algorithms, which can be specified by passing the named argument `strategy=...`.

| Strategy Name | Notes | Interface |
|---|---|---|
| `RecoStrategy.Best` | Dynamically switch strategy based on input particle density | `jet_reconstruct` |
| `RecoStrategy.N2Plain` | Global matching of particles at each interation (works well for low $N$) | `plain_jet_reconstruct` |
| `RecoStrategy.N2Tiled` | Use tiles of radius $R$ to limit search space (works well for higher $N$) | `tiled_jet_reconstruct` |

Generally one can use the `jet_reconstruct` interface, shown above, as the
*Best* strategy safely as the overhead is extremely low. That interface supports
a `strategy` option to switch to a different option.

For ``e^+e^-`` algorithms particle densities are low, so the only implementation
is of the same type as `N2Plain`.

## Inclusive and Exclusive Selections

To obtain final jets both inclusive (``p_T`` cut) and exclusive (``n_{jets}`` or
``d_{ij}`` cut) selections are supported:

- [inclusive_jets(clusterseq::ClusterSequence, ptmin = 0.0)](@ref)
- [exclusive_jets(clusterseq::ClusterSequence; dcut = nothing, njets = nothing)](@ref)

### Sorting

Sorting vectors is trivial in Julia, no special sorting methods are provided. As
an example, to sort exclusive jets of ``>5.0`` (usually GeV, depending on your
EDM) from highest energy to lowest:

```julia
sorted_jets = sort!(inclusive_jets(cs::ClusterSequence; ptmin=5.0), 
  by=JetReconstruction.energy, rev=true)
```

## Plotting and Animation

![illustration](assets/jetvis.png)

To visualise the clustered jets as a 3d bar plot (see illustration above) we now
use `Makie.jl`. See the `jetsplot` function in `ext/JetVisualisation.jl` and its
documentation for more. There are two worked examples in the `examples`
directory.

The plotting code is a package extension and will load if the one of the `Makie`
modules is loaded in the environment.

The [`animatereco`](@ref) function will animate the reconstruction sequence, given a
`ClusterSequence` object. See the function documentation for the many options
that can be customised.

## Serialisation

The package also provides methods such as `loadjets`, `loadjets!`, and
`savejets` that one can use to save and load objects on/from disk easily in a
very flexible format. See documentation for more.

## Reference

Although it has been developed further since the CHEP2023 conference, the CHEP
conference proceedings,
[10.1051/epjconf/202429505017](https://doi.org/10.1051/epjconf/202429505017),
should be cited if you use this package:

```bibtex
@article{refId0,
    author = {{Stewart, Graeme Andrew} and {Gras, Philippe} and {Hegner, Benedikt} and {Krasnopolski, Atell}},
    doi = {10.1051/epjconf/202429505017},
    journal = {EPJ Web of Conf.},
    pages = {05017},
    title = {Polyglot Jet Finding},
    url = {https://doi.org/10.1051/epjconf/202429505017},
    volume = 295,
    year = 2024,
    eprint={2309.17309},
    archivePrefix={arXiv},
    primaryClass={hep-ex}
}
```

The original paper on [arXiv](https://arxiv.org/abs/2309.17309) is:

```bibtex
@misc{stewart2023polyglot,
      title={Polyglot Jet Finding}, 
      author={Graeme Andrew Stewart and Philippe Gras and Benedikt Hegner and Atell Krasnopolski},
      year={2023},
      eprint={2309.17309},
      archivePrefix={arXiv},
      primaryClass={hep-ex}
}
```

## Authors and Copyright

Code in this package is authored by:

- Atell Krasnopolski <delta_atell@protonmail.com>
- Graeme A Stewart <graeme.andrew.stewart@cern.ch>
- Philippe Gras <philippe.gras@cern.ch>

and is Copyright 2022-2024 The Authors, CERN.

The code is under the MIT License.
