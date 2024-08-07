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
which are used in high-energy physics as part of event reconstruction for $pp$
and $e^+e^-$ colliders.

## Algorithms

Algorithms used are based on the C++ FastJet package (<https://fastjet.fr>,
[hep-ph/0512210](https://arxiv.org/abs/hep-ph/0512210),
[arXiv:1111.6097](https://arxiv.org/abs/1111.6097)), reimplemented natively in
Julia.

The algorithms include ``\text{anti}-{k}_\text{T}``, Cambridge/Aachen and
inclusive ``k_\text{T}``.

## Reconstruction Interface

The main interface for reconstruction is [`jet_reconstruct`](@ref), called as, e.g.,

```julia
jet_reconstruct(particles; p = -1, R = 1.0)
```

or

```julia
jet_reconstruct(particles; algorithm = JetAlgorithm.AntiKt, R = 1.0)
```

Where `particles` is a collection of 4-vector objects to reconstruct and the
algorithm is either given explicitly or implied by the power value. For the case of generalised $k_T$ both the algorithm (`JetAlgorithm.GenKt`) and `p` are needed.

The object returned is a [`ClusterSequence`](@ref), which internally tracks all
merge steps.

## Strategy

Three strategies are available for the different algorithms:

| Strategy Name | Notes | Interface |
|---|---|---|
| `RecoStrategy.Best` | Dynamically switch strategy based on input particle density | `jet_reconstruct` |
| `RecoStrategy.N2Plain` | Global matching of particles at each interation (works well for low $N$) | `plain_jet_reconstruct` |
| `RecoStrategy.N2Tiled` | Use tiles of radius $R$ to limit search space (works well for higher $N$) | `tiled_jet_reconstruct` |

Generally one can use the `jet_reconstruct` interface, shown above, as the
*Best* strategy safely as the overhead is extremely low. That interface supports
a `strategy` option to switch to a different option.

Another option, if one wishes to use a specific strategy, is to call that
strategy's interface directly, e.g.,

```julia
# For N2Plain strategy called directly
plain_jet_reconstruct(particles::Vector{T}; p = -1, R = 1.0, recombine = +)
```

Note that there is no `strategy` option in these interfaces.

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
