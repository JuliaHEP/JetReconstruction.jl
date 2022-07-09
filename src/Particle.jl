"""
This module defines the basic approach to representing particles and pseudojets. That is, in the form of an indexable (either mutable or immutable) object `obj` where
```
obj[1] # energy
obj[2] # px
obj[3] # py
obj[4] # pz
```
Therefore, a simple vector [E, x, y, z] or a [static array](https://github.com/JuliaArrays/StaticArrays.jl) will do fine.

You can define your own specific data type to represent physical objects. For the anti-kt algorithm to work with a data type, for instance, it only needs to have `pt`, `phi`, `eta`, and `+` operations defined.
"""
module Particle

# TODO: add built-in integration with LorentzVector and LorentzVectorHEP and write about it in the docstring

export energy, px, py, pz, pt, phi, mass, eta, kt, ϕ, η

@inline energy(p) = p[1]

@inline px(p) = p[2]

@inline py(p) = p[3]

@inline pz(p) = p[4]

@inline pt(p) = @fastmath sqrt(p[2]^2 + p[3]^2)
kt = pt

@inline phi(p) = @fastmath atan(p[3], p[2])
ϕ = phi

@inline mass(p) = @fastmath sqrt(max(p[0]^2 - p[2]^2 - p[3]^2 - p[4]^2, 0))

@inline eta(p) = asinh(p[4]/pt(p))
η = eta

end
