"""
This module defines the basic approach to particle and pseudojet representation. That is, in the form of an indexable (either mutable or immutable) object `obj` where
```
obj[1] # px
obj[2] # py
obj[3] # pz
obj[4] # energy
```
Therefore, a simple vector [x, y, z, E] or a [static array](https://github.com/JuliaArrays/StaticArrays.jl) will do fine.

You can define your own specific data type to represent physical objects. For the anti-kt algorithm to work with a data type, for instance, it only needs to have `pt`, `phi`, `eta`, and `+` operations defined. Be sure, however, to define those functions as extentions to the already existing `JetReconstruction.pt`, `JetReconstruction.phi`, `JetReconstruction.eta` and `Base.:+` respectively:
```
import JetReconstruction

JetReconstruction.eta(x::YourType) = # your eta definition
JetReconstruction.phi(x::YourType) = # your phi definition
JetReconstruction.pt(x::YourType) = # your pt definition

Base.:+(x::YourType, y::YourType) = # your + definition (only used as a recombination function)
```
Since `+` is only used in recombination, you can leave it undefined, if you use a custom recombination routine.
"""

"""
"""
@inline energy(p) = p[4]

@inline px(p) = p[1]

@inline py(p) = p[2]

@inline pz(p) = p[3]

@inline pt2(p) = @fastmath px(p)^2 + py(p)^2
@inline pt(p) = @fastmath sqrt(pt2(p))
const kt = pt

# @inline phi(p) = @fastmath atan(py(p), px(p))
# Fix to return in range [0, 2π)
phi(p) = begin
    if pt2(p) == 0.0
        phi = 0.0
    else
        phi = atan(py(p), px(p))
    end
    if phi < 0.0
        phi += 2π
    end
    phi
end
const ϕ = phi

@inline mass(p) = @fastmath sqrt(max(energy(p)^2 - px(p)^2 - py(p)^2 - pz(p)^2, 0))

#@inline pseudorap(p) = asinh(pz(p)/pt(p)) # pseudorapidity

# Rapidity
eta(p) = begin
    _pt2 = pt2(p)
    _abspz = abs(pz(p))
    if (energy(p) == _abspz) && (_pt2 == 0.0)
        return (-1)^(pz(p) < 0)*(1e5 + _abspz) # a very large value that depends on pz
    end
    _m2 = max((energy(p) + pz(p))*(energy(p) - pz(p)) - _pt2, 0.0) # mass^2
    E_plus_z = energy(p) + _abspz
    return (-1)^(pz(p) > 0) * 0.5*log((_pt2 + _m2)/(E_plus_z^2))
end
const η = eta
