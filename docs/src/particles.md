# Input Particle Types

For the `particles` input to the reconstruction any one dimensional
`AbstractArray{T, 1}` can be used, where the type `T` has to implement methods
to extract the 4-vector components, viz, the following are required:

- `JetReconstuction.px(particle::T)`
- `JetReconstuction.py(particle::T)`
- `JetReconstuction.pz(particle::T)`
- `JetReconstuction.energy(particle::T)`

Currently built-in supported types are
[`LorentzVectorHEP`](https://github.com/JuliaHEP/LorentzVectorHEP.jl), the
`PseudoJet` and `EEJet`s from this package, and `ReconstructedParticles` from
[EDM4hep Inputs](@ref).

If you require support for a different input collection type then ensure you
define the `px()`, etc. methods *for your specific type* and *in the
`JetReconstruction` package*. This use of what might be considered type piracy
is blessed as long as you are en *end user* of the jet reconstruction package.

If your type is used in several places or by different users, please consider
writing a package extension that will support your type, following the model for
EDM4hep in `ext/EDM4hepJets.jl`.
