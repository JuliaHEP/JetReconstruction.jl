# Input Particle Types

For the `particles` input to the reconstruction any one dimensional
`AbstractArray{T, 1}` can be used, in general the type `T` should implement the
[LorentzVectorBase](https://github.com/JuliaHEP/LorentzVectorBase.jl)
interface. That is the `T` should be a valid coordinate system and the
following methods are required:

- `LorentzVectorBase.px(particle::T)`
- `LorentzVectorBase.py(particle::T)`
- `LorentzVectorBase.pz(particle::T)`
- `LorentzVectorBase.energy(particle::T)`

In addition, direct support for the construction of particles from
`LorentzVector` and `LorentzVectorCyl` types is provided.

Direct use of the `PseudoJet` and `EEJet`s from this package, and
`ReconstructedParticles` from [EDM4hep Inputs](@ref) is supported.
