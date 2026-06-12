# Inspired by the PseudoJet class of c++ code of Fastjet (https://fastjet.fr,
#  hep-ph/0512210,  arXiv:1111.6097)
#
# Some of the implementation is taken from LorentzVectorHEP.jl, by Jerry Ling

import Accessors

"""
    struct PseudoJet{T <: Real} <: FourMomentum{T}

The `PseudoJet` struct represents a pseudojet, a four-momentum object used in
jet reconstruction algorithms. It is parameterized by the type `T` of its 
momentum components.

Additional information for the link back into the history of the clustering is 
stored in the `_cluster_hist_index` field. There is caching of the more 
expensive calculations for rapidity and azimuthal angle.

# Fields
- `px::T`: The x-component of the momentum.
- `py::T`: The y-component of the momentum.
- `pz::T`: The z-component of the momentum.
- `E::T`: The energy component of the momentum.
- `_cluster_hist_index::Int`: The index of the cluster history that corresponds 
  to this PseudoJet
- `_pt2::T`: The squared transverse momentum.
- `_rap::T`: The rapidity.
- `_phi::T`: The azimuthal angle.

"""
struct PseudoJet{T <: Real} <: FourMomentum{T}
    px::T
    py::T
    pz::T
    E::T
    _cluster_hist_index::Int
    _pt2::T
    _rap::T
    _phi::T
end

"""
    PseudoJet(px::T, py::T, pz::T, E::T; cluster_hist_index::Integer = 0) where {T <: Real}

Construct a parameterised `PseudoJet{T}` from a four momentum `(px, py, pz, E)`
with cluster index `cluster_hist_index`.

# Details

If the (default) value of `cluster_hist_index=0` is used, the PseudoJet cannot be
used in a reconstruction sequence.
"""
function PseudoJet(px::T, py::T, pz::T, E::T;
                   cluster_hist_index::Integer = 0) where {T <: Real}
    @muladd pt2 = px * px + py * py
    phi = pt2 == 0.0 ? 0.0 : atan(py, px)
    phi = phi < 0.0 ? phi + 2π : phi
    if E == abs(pz) && iszero(pt2)
        # Point has infinite rapidity - convert that into a very large
        #    number, but in such a way that different 0-pt momenta will have
        #    different rapidities (so as to lift the degeneracy between
        #    them) [this can be relevant at parton-level]
        MaxRapHere = _MaxRap + abs(pz)
        rap = pz >= 0.0 ? MaxRapHere : -MaxRapHere
    else
        # get the rapidity in a way that's modestly insensitive to roundoff
        # error when things pz,E are large (actually the best we can do without
        # explicit knowledge of mass)
        effective_m2 = (E + pz) * (E - pz) - pt2
        effective_m2 = max(0.0, effective_m2) # force non tachyonic mass
        E_plus_pz = E + abs(pz) # the safer of p+, p-
        rap = 0.5 * log((pt2 + effective_m2) / (E_plus_pz * E_plus_pz))
        rap = pz > 0.0 ? -rap : rap
    end
    PseudoJet{T}(px, py, pz, E, cluster_hist_index, pt2, rap, phi)
end

"""
    PseudoJet(px, py, pz, E; cluster_hist_index::Int = 0)

Constructor for mixed input numerical types. The inputs are promoted to a
common type `T`, and a `PseudoJet{T}` is returned.
"""
function PseudoJet(px::Tpx, py::Tpy, pz::Tpz, E::TE;
                   cluster_hist_index::Int = 0) where {Tpx, Tpy, Tpz, TE}
    PseudoJet(promote(px, py, pz, E)...; cluster_hist_index)
end

"""
    PseudoJet{T}(px, py, pz, E; cluster_hist_index::Integer = 0) where {T <: Real}

Constructor for explicit parameter type `T`.
"""
function PseudoJet{T}(px, py, pz, E; cluster_hist_index::Integer = 0) where {T <: Real}
    PseudoJet(T(px), T(py), T(pz), T(E); cluster_hist_index)
end

"""
    const invalid_pseudojet::PseudoJet{Float64}

Used to mark an invalid result in case the corresponding substructure tagging fails.
This constant uses `Float64` as the default type.

## Deprecated
This value will be removed in future. It should not be used directly, instead use the `isvalid()` function.
"""
const invalid_pseudojet = PseudoJet(0.0, 0.0, 0.0, 0.0; cluster_hist_index = typemin(Int))

import Base.zero
"""
    zero(::Type{PseudoJet{T}}) where {T <: Real}

Generate an invalid `PseudoJet{T}`, used to mark an invalid result in case the 
corresponding substructure tagging fails.

## Note
For technical reasons the rapidity is set to `_MaxRap`.
"""
function zero(::Type{PseudoJet{T}}) where {T <: Real}
    PseudoJet{T}(0.0, 0.0, 0.0, 0.0, typemin(Int), 0.0, _MaxRap, 0.0)
end

"""
    PseudoJet{T}(;pt, rap, phi, m = 0, cluster_hist_index::Int = 0) where {T <: Real}

Construct a `PseudoJet{T}` from `(pt, rap, phi, m)` with the cluster index
`cluster_hist_index`.

# Details

If the (default) value of `cluster_hist_index=0` is used, the PseudoJet cannot be
used in a reconstruction sequence.
"""
function PseudoJet{T}(; pt, rap, phi, m = 0,
                      cluster_hist_index::Int = 0) where {T <: Real}
    phi = phi < 0 ? phi + 2π : phi
    phi = phi > 2π ? phi - 2π : phi
    ptm = (m == 0) ? pt : sqrt(pt^2 + m^2)
    exprap = exp(rap)
    pminus = ptm / exprap
    pplus = ptm * exprap
    px = pt * cos(phi)
    py = pt * sin(phi)
    pz = @fastmath (pplus - pminus) / 2
    E = @fastmath (pplus + pminus) / 2

    (px, py, pz, E, pt, rap, phi) = promote(px, py, pz, E, pt, rap, phi)

    PseudoJet{T}(px, py, pz, E, cluster_hist_index, pt^2, rap, phi)
end

function PseudoJet(; pt::T, rap::T, phi::T, m = 0,
                   cluster_hist_index::Int = 0) where {T <: Real}
    PseudoJet{T}(pt = pt, rap = rap, phi = phi, m = m,
                 cluster_hist_index = cluster_hist_index)
end

"""
    PseudoJet{T}(jet::PseudoJet{T}; cluster_hist_index::Int = 0) where {T <: Real}

Construct a `PseudoJet` from another `PseudoJet` object and assign the given cluster index to it.
"""
function PseudoJet{T}(jet::PseudoJet{T}; cluster_hist_index::Int = 0) where {T <: Real}
    Accessors.@set jet._cluster_hist_index = cluster_hist_index
end

"""
    PseudoJet{O}(jet::PseudoJet{T}; cluster_hist_index::Int = 0) where {T <: Real, O <: Real}

Construct a `PseudoJet{O}` from another `PseudoJet{T}` object and assign given
cluster index to it. Used to change the underlying type from `T` to `O`.
"""
function PseudoJet{O}(jet::PseudoJet{T};
                      cluster_hist_index::Int = 0) where {T <: Real, O <: Real}
    PseudoJet{O}(jet.px, jet.py, jet.pz, jet.E; cluster_hist_index)
end

"""
    PseudoJet(jet::LorentzVector; cluster_hist_index::Int = 0)

Construct a `PseudoJet` from a `LorentzVector` object with the cluster index.
The underlying type is inferred from the `LorentzVector` components.
"""
function PseudoJet(jet::LorentzVector; cluster_hist_index::Int = 0)
    PseudoJet(jet.x, jet.y, jet.z, jet.t; cluster_hist_index)
end

"""
    PseudoJet{T}(jet::LorentzVector; cluster_hist_index::Int = 0) where {T <: Real}

Construct a `T` parameterised `PseudoJet` from a `LorentzVector` object with
the cluster index.
"""
function PseudoJet{T}(jet::LorentzVector; cluster_hist_index::Int = 0) where {T <: Real}
    PseudoJet(T(jet.x), T(jet.y), T(jet.z), T(jet.t); cluster_hist_index)
end

"""
    PseudoJet(jet::LorentzVectorCyl; cluster_hist_index::Int = 0)

Construct a `PseudoJet` from a `LorentzVectorCyl` object with the given cluster index.
The underlying type is inferred from the `LorentzVectorCyl` components.
"""
function PseudoJet(jet::LorentzVectorCyl; cluster_hist_index::Int = 0)
    PseudoJet(; pt = pt(jet), rap = rapidity(jet), phi = phi(jet), m = mass(jet),
              cluster_hist_index)
end

"""
function PseudoJet{T}(jet::LorentzVectorCyl; cluster_hist_index::Int = 0) where {T <: Real}

Construct a `T` parameterised `PseudoJet` from a `LorentzVectorCyl` object with
the given cluster index.
"""
function PseudoJet{T}(jet::LorentzVectorCyl; cluster_hist_index::Int = 0) where {T <: Real}
    PseudoJet(; pt = T(pt(jet)), rap = T(rapidity(jet)), phi = T(phi(jet)),
              m = T(mass(jet)),
              cluster_hist_index)
end

"""
    PseudoJet(jet::Any; cluster_hist_index::Int = 0)

Construct a `PseudoJet` from a generic object `jet` with the given cluster index.
The generic jet must implement the `LorentzVectorBase` interface.

# Details

This function is used to convert a generic object `jet` into a `PseudoJet`.
The generic jet should implement the `LorentzVectorBase` interface: the
`LorentzVectorBase.{px,py,pz,energy}()` methods will be called to retrieve
the four momentum components.

The `cluster_hist_index` is optional, but needed if the `jet` is part of a
reconstruction sequence. If not provided, it defaults to `0` as an "invalid"
value.
"""
function PseudoJet(jet::Any; cluster_hist_index::Int = 0)
    # Check that the interface is implemented
    if hasmethod(LorentzVectorBase.coordinate_system, (typeof(jet),))
        return PseudoJet(LorentzVectorBase.px(jet), LorentzVectorBase.py(jet),
                         LorentzVectorBase.pz(jet), LorentzVectorBase.energy(jet);
                         cluster_hist_index)
    else
        throw(ArgumentError("PseudoJet cannot be constructed from object of type '$(typeof(jet))'"))
    end
end

"""
    PseudoJet{T}(jet::Any; cluster_hist_index::Int = 0)

Construct a `PseudoJet{T}` from a generic object `jet` with the given cluster
index. The generic jet must implement the `LorentzVectorBase` interface.

# Details

This function is used to convert a generic object `jet` into a `PseudoJet`.
The generic jet should implement the `LorentzVectorBase` interface: the
`LorentzVectorBase.{px,py,pz,energy}()` methods will be called to retrieve
the four momentum components.

The `cluster_hist_index` is optional, but needed if the `jet` is part of a
reconstruction sequence. If not provided, it defaults to `0` as an "invalid"
value.
"""
function PseudoJet{T}(jet::Any; cluster_hist_index::Int = 0) where {T <: Real}
    # Check that the interface is implemented
    if hasmethod(LorentzVectorBase.coordinate_system, (typeof(jet),))
        return PseudoJet(T(LorentzVectorBase.px(jet)), T(LorentzVectorBase.py(jet)),
                         T(LorentzVectorBase.pz(jet)), T(LorentzVectorBase.energy(jet));
                         cluster_hist_index)
    else
        throw(ArgumentError("PseudoJet cannot be constructed from object of type '$(typeof(jet))'"))
    end
end

import Base.isvalid
"""
    isvalid(j::PseudoJet)

Function to check whether a given `PseudoJet` object is non-zero or not.
Primarily used for checking the validity of outputs of substructure tagging.

# Returns
- `Bool`: `true` if the `PseudoJet` object is non-zero (valid), `false` otherwise. 
"""
@inline isvalid(j::PseudoJet) = !(cluster_hist_index(j) == typemin(Int))

"""
    phi(p::PseudoJet)

Return the azimuthal angle, ϕ, of a `PseudoJet` object `p` in the range
[0, 2π). This accessor uses the pre-calculated value that the struct has.

Note that the range [0, 2π) differs from the convention in `LorentzVectorBase`,
which is [-π, π].
"""
phi(p::PseudoJet) = p._phi

"""
    rapidity(p::PseudoJet)

Return the rapidity of a `PseudoJet` object. This accessor uses the
pre-calculated value that the struct has.
"""
rapidity(p::PseudoJet) = p._rap
LorentzVectorBase.rapidity(p::PseudoJet) = p._rap

"""
    pt2(p::PseudoJet)

Return the squared transverse momentum of a `PseudoJet`. This accessor uses the
pre-calculated value that the struct has.
"""
pt2(p::PseudoJet) = p._pt2
LorentzVectorBase.pt2(p::PseudoJet) = p._pt2

"""
    pt(p::PseudoJet)

Return the scalar transverse momentum (pt) of a `PseudoJet`. This accessor uses
the pre-calculated value that the struct has.
"""
pt(p::PseudoJet) = sqrt(p._pt2)
LorentzVectorBase.pt(p::PseudoJet) = sqrt(p._pt2)
