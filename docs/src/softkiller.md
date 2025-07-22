# SoftKiller

The **SoftKiller** algorithm is a pileup mitigation technique for jet reconstruction in high-energy physics. It is designed to suppress the effect pileup by removing low-``p_T`` particles before jet clustering.

## Overview

SoftKiller divides the rapidity-azimuth (``y\text{-}\phi``) plane into a regular grid. In each grid cell, it determines a ``p_T`` threshold such that at least one particle per cell survives. All particles with ``p_T`` below this threshold are removed. This approach efficiently reduces pileup contamination while preserving the hard event structure.

The method is described in [Cacciari, Salam, Soyez, Eur. Phys. J. C75 (2015) 59](https://arxiv.org/abs/1407.0408).

## Interface

The main interface for SoftKiller is the [`softkiller`](@ref) function, which is used as follows:

```julia
using JetReconstruction

sk = SoftKiller(rapmax=5.0, grid_size=0.4)
filtered_particles, pt_threshold = softkiller(sk, particles)
```

- `particles` is a collection of `PseudoJet` objects (see [Input Particle Types](@ref)).
- `rapmax` sets the maximum rapidity considered.
- `grid_size` sets the size of each grid cell in rapidity and azimuthal angle.

The function returns the filtered list of particles and the computed ``p_T`` threshold.

## Example

A typical workflow using SoftKiller is:

```julia
using JetReconstruction

# Read event particles (see examples/softkiller/softkiller_plots.jl)
particles = read_final_state_particles("event.hepmc3", PseudoJet)[1]

# Set up SoftKiller
sk = SoftKiller(rapmax=5.0, grid_size=0.4)

# Apply SoftKiller
filtered_particles, pt_threshold = softkiller(sk, particles)

# Cluster jets as usual
cs = jet_reconstruct(filtered_particles; algorithm=JetAlgorithm.AntiKt, R=0.4)
jets = inclusive_jets(cs, ptmin=25.0)
```

See [`examples/softkiller/softkiller_plots.jl`] and [`examples/softkiller/softkiller_runtime.jl`] for more specific working examples.

## Parameters

| Parameter   | Description                                   | Default |
|-------------|-----------------------------------------------|---------|
| `rapmax`    | Maximum rapidity for grid                     | 5.0     |
| `grid_size` | Size of grid cell in rapidity-azimuth         | 0.4     |

## Integration

SoftKiller is typically used as a preprocessing step before jet clustering. It is compatible with all jet algorithms provided by this package and is especially useful in high pileup environments.

## References

- M. Cacciari, G. P. Salam, G. Soyez, [Eur. Phys. J. C75 (2015) 59](https://arxiv.org/abs/1407.0408)
- [SoftKiller in FastJet](https://fastjet.fr/contrib/) and (https://phab.hepforge.org/source/fastjetsvn/browse/contrib/contribs/SoftKiller/tags/1.0.0/)

## Authors

The SoftKiller implementation in this package follows the approach described in the original paper and is integrated with the native Julia jet reconstruction algorithms.
