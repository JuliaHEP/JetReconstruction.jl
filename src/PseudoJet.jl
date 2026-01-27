# Inspired by the PseudoJet class of c++ code of Fastjet (https://fastjet.fr,
#  hep-ph/0512210,  arXiv:1111.6097)
#
# Some of the implementation is taken from LorentzVectorHEP.jl, by Jerry Ling

import Accessors

"""
    struct PseudoJet <: FourMomentum

The `PseudoJet` struct represents a pseudojet, a four-momentum object used in
jet reconstruction algorithms. Additional information for the link back into the
history of the clustering is stored in the `_cluster_hist_index` field. There is
caching of the more expensive calculations for rapidity and azimuthal angle.

# Fields
- `px::Float64`: The x-component of the momentum.
- `py::Float64`: The y-component of the momentum.
- `pz::Float64`: The z-component of the momentum.
- `E::Float64`: The energy component of the momentum.
- `_cluster_hist_index::Int`: The index of the cluster history.
- `_pt2::Float64`: The squared transverse momentum.
- `_rap::Float64`: The rapidity.
- `_phi::Float64`: The azimuthal angle.

"""
struct PseudoJet <: FourMomentum
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
    _cluster_hist_index::Int
    _pt2::Float64
    _rap::Float64
    _phi::Float64
end

"""
    PseudoJet(px::Real, py::Real, pz::Real, E::Real; cluster_hist_index::Int = 0)

Construct a PseudoJet from a four momentum `(px, py, pz, E)`` with cluster index
`cluster_hist_index`.

# Details

If the (default) value of `cluster_hist_index=0` is used, the PseudoJet cannot be
used in a reconstruction sequence.
"""
function PseudoJet(px::Real, py::Real, pz::Real, E::Real; cluster_hist_index::Int = 0)
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
    PseudoJet(px, py, pz, E, cluster_hist_index, pt2, rap, phi)
end

"""
    const invalid_pseudojet = PseudoJet(0.0, 0.0, 0.0, 0.0)

Used to mark an invalid result in case the corresponding substructure tagging fails."""
const invalid_pseudojet = PseudoJet(0.0, 0.0, 0.0, 0.0)

"""
    PseudoJet(;pt::Real, rap::Real, phi::Real, m::Real = 0, cluster_hist_index::Int = 0)

Construct a PseudoJet from `(pt, y, ϕ, m)` with the cluster index
`cluster_hist_index`.

# Details

If the (default) value of `cluster_hist_index=0` is used, the PseudoJet cannot be
used in a reconstruction sequence.
"""
function PseudoJet(; pt::Real, rap::Real, phi::Real, m::Real = 0,
                   cluster_hist_index::Int = 0)
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

    PseudoJet(px, py, pz, E, cluster_hist_index, pt^2, rap, phi)
end

"""
    PseudoJet(jet::PseudoJet; cluster_hist_index::Int = 0)

Construct a PseudoJet from another `PseudoJet` object and assign given cluster index to it.
"""
function PseudoJet(jet::PseudoJet; cluster_hist_index::Int = 0)
    Accessors.@set jet._cluster_hist_index = cluster_hist_index
end

"""
    PseudoJet(jet::LorentzVector; cluster_hist_index::Int = 0)

Construct a PseudoJet from a `LorentzVector` object with the cluster index.
"""
function PseudoJet(jet::LorentzVector; cluster_hist_index::Int = 0)
    PseudoJet(jet.x, jet.y, jet.z, jet.t; cluster_hist_index)
end

"""
    PseudoJet(jet::LorentzVectorCyl; cluster_hist_index::Int = 0)

Construct a PseudoJet from a `LorentzVectorCyl` object with the given cluster index.
"""
function PseudoJet(jet::LorentzVectorCyl; cluster_hist_index::Int = 0)
    PseudoJet(; pt = pt(jet), rap = rapidity(jet), phi = phi(jet), m = mass(jet),
              cluster_hist_index)
end

"""
    PseudoJet(jet::Any; cluster_hist_index::Int = 0)

Construct a PseudoJet from a generic object `jet` with the given cluster index.
These generic jets must implement the LorentzVectorBase interface from
`LorentzVectorBase.jl`.

# Details

This function is used to convert a generic object `jet` into a `PseudoJet`.
These generic jets should implement the LorentzVectorBase interface: the
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

import Base.isvalid
"""
    isvalid(j::PseudoJet)

Function to check whether a given `PseudoJet` object is non-zero or not.
Primarily to use for checking the validity of outputs of substructure tagging.

# Returns
- `Bool`: `true` if the `PseudoJet` object is non-zero (valid), `false` otherwise. 
"""
isvalid(j::PseudoJet) = !(j === invalid_pseudojet)

"""
    phi(p::PseudoJet)

Return the azimuthal angle, ϕ, of a `PseudoJet` object `p` in the range
[0, 2π). This accessor uses the pre-calculated value that the struct has.
"""
phi(p::PseudoJet) = p._phi

"""
    rapidity(p::PseudoJet)

Return the rapidity of a `PseudoJet` object. This accessor uses the
pre-calculated value that the struct has.
"""
rapidity(p::PseudoJet) = p._rap

"""
    pt2(p::PseudoJet)

Return the squared transverse momentum of a `PseudoJet`. This accessor uses the
pre-calculated value that the struct has.
"""
pt2(p::PseudoJet) = p._pt2

"""
    pt(p::PseudoJet)

Return the scalar transverse momentum (pt) of a PseudoJet. This accessor uses
the precalculated value that the struct has.
"""
pt(p::PseudoJet) = sqrt(p._pt2)
