# JetReconstruction.jl

[![Build Status](https://github.com/JuliaHEP/JetReconstruction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaHEP/JetReconstruction.jl/actions/workflows/CI.yml?query=branch%3Amain)

## This package implements sequential Jet Reconstruction (clustering)

### Algorithms

Algorithms used are based on the C++ FastJet package (<https://fastjet.fr>,
[hep-ph/0512210](https://arxiv.org/abs/hep-ph/0512210),
[arXiv:1111.6097](https://arxiv.org/abs/1111.6097)), reimplemented natively in Julia.

One algorithm is the plain N2 approach `N2Plain`, the other uses the FastJet tiling
approach, `N2Tiled`.

### Interface

To use the package one should call the appropriate API interface of the desired algorithm

```julia
# N2Plain
plain_jet_reconstruct(particles::Vector{T}; p = -1, R = 1.0, recombine = +, ptmin = 0.0)
```

```julia
# N2Tiled
tiled_jet_reconstruct(particles::Vector{T}; p = -1, R = 1.0, recombine = +, ptmin = 0.0)
```

- `particles` - a vector of input particles for the clustering
  - Any type that supplies the methods `pt2()`, `phi()`, `rapidity()`, `px()`, `py()`, `pz()`, `energy()` can be used
  - These methods have to be defined in the namespace of this package, i.e., `JetReconstruction.pt2(::T)`
  - The `PseudoJet` type from this package, or a 4-vector from `LorentzVectorHEP` are suitable (and have the appropriate definitions)
- `p` - the transverse momentum power used in the $d_{ij}$ metric for deciding on closest jets, as $k^{2p}_\text{T}$
  - `-1` gives anti-${k}_\text{T}$ clustering
  - `0` gives Cambridge/Achen
  - `1` gives inclusive $k_\text{T}$
- `R` - the cone size parameter (no particles more geometrically distance than `R` will be merged)
- `recombine` - the function used to merge two pseudojets (by default this is simple 4-vector addition of $(E, \mathbf{p})$)
- `ptmin` - only jets of transverse momentum greater than or equal to `ptmin` will be returned

Both of these functions return tuple of `(jets, seq)`, where these are a vector of `JetReconstruction.PseudoJet` objects and cluster sequencing information, respectively.

Note that the `N2Plain` algorithm is usually faster at low densities, $\lessapprox 50$ input particles, otherwise the tiled algorithm is faster.

### Example

See the `examples/jetreco.jl` script for a full example of how to call jet reconstruction.

```julia
julia --project=. examples/jetreco.jl --maxevents=100 --nsamples=1 --strategy=N2Plain test/data/events.hepmc3
...
julia --project=. examples/jetreco.jl --maxevents=100 --nsamples=1 --strategy=N2Tiled test/data/events.hepmc3
...
```

The example also shows how to use `JetReconstruction.HepMC3` to read HepMC3 ASCII files (via the `read_final_state_particles()` wrapper).

### Plotting

**TO BE FIXED**

![illustration](img/illustration.jpeg)

To visualise the clustered jets as a 3d bar plot (see illustration above) we now use `Makie.jl`. See the `jetsplot` function and its documentation for more. 

## Reference

Although it has been developed further since the CHEP2023 conference, the CHEP conference proceedings, [arXiv:2309.17309](https://arxiv.org/abs/2309.17309), should be cited if you use this package:

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

and is Copyright 2022-2023 The Authors, CERN.

The code is under the MIT License.
