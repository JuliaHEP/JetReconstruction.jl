# Adapted from PseudoJet class of c++ code of Fastjet (https://fastjet.fr,
#  hep-ph/0512210,  arXiv:1111.6097)
#
#  Copyright (c) 2005-2020, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
#
# Some of the implementation is taken from LorentzVectorHEP.jl, (c) Jerry Ling

"""Interface for composite types that includes fields px, py, py, and E
that represents the components of a four-momentum vector."""
abstract type FourMomentum end

"""Used to protect against parton-level events where pt can be zero
for some partons, giving rapidity=infinity. KtJet fails in those cases."""
const _MaxRap = 1e5

"""Used to signal that the phi value has not yet been computed."""
const _invalid_phi = -100.0

"""Used to signal that the rapidity value has not yet been computed."""
const _invalid_rap = -1.e200

# @ingroup basic_classes
# \class PseudoJet
# Class to contain pseudojets, including minimal information of use to
# jet-clustering routines.

"""
    mutable struct PseudoJet <: FourMomentum

The `PseudoJet` struct represents a pseudojet, a four-momentum object used in
jet reconstruction algorithms. Additonal information for the link back into the
history of the clustering is stored in the `_cluster_hist_index` field. There is
caching of the more expensive calculations for rapidity and azimuthal angle.

# Fields
- `px::Float64`: The x-component of the momentum.
- `py::Float64`: The y-component of the momentum.
- `pz::Float64`: The z-component of the momentum.
- `E::Float64`: The energy component of the momentum.
- `_cluster_hist_index::Int`: The index of the cluster history.
- `_pt2::Float64`: The squared transverse momentum.
- `_inv_pt2::Float64`: The inverse squared transverse momentum.
- `_rap::Float64`: The rapidity.
- `_phi::Float64`: The azimuthal angle.

"""
mutable struct PseudoJet <: FourMomentum
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
    _cluster_hist_index::Int
    _pt2::Float64
    _inv_pt2::Float64
    _rap::Float64
    _phi::Float64
end

"""
    PseudoJet(px::Float64, py::Float64, pz::Float64, E::Float64,
        _cluster_hist_index::Int,
        pt2::Float64)

Constructs a PseudoJet object with the given momentum components and energy and
history index.

# Arguments
- `px::Float64`: The x-component of the momentum.
- `py::Float64`: The y-component of the momentum.
- `pz::Float64`: The z-component of the momentum.
- `E::Float64`: The energy.
- `_cluster_hist_index::Int`: The cluster history index.
- `pt2::Float64`: The transverse momentum squared.

# Returns
A `PseudoJet` object.
"""
PseudoJet(px::Float64, py::Float64, pz::Float64, E::Float64,
_cluster_hist_index::Int,
pt2::Float64) = PseudoJet(px,
                          py, pz, E, _cluster_hist_index,
                          pt2, 1.0 / pt2, _invalid_rap, _invalid_phi)

"""
    PseudoJet(px::Float64, py::Float64, pz::Float64, E::Float64)

Constructs a PseudoJet object with the given momentum components and energy.

# Arguments
- `px::Float64`: The x-component of the momentum.
- `py::Float64`: The y-component of the momentum.
- `pz::Float64`: The z-component of the momentum.
- `E::Float64`: The energy.

# Returns
A PseudoJet object.
"""
PseudoJet(px::Float64, py::Float64,
pz::Float64, E::Float64) = PseudoJet(px, py, pz, E, 0, px^2 + py^2)

import Base.show
"""
    show(io::IO, jet::PseudoJet)

Print a `PseudoJet` object to the specified IO stream.

# Arguments
- `io::IO`: The IO stream to which the information will be printed.
- `jet::PseudoJet`: The `PseudoJet` object whose information will be printed.
"""
show(io::IO, jet::PseudoJet) = begin
    print(io, "Pseudojet(px: ", jet.px, " py: ", jet.py, " pz: ", jet.pz, " E: ", jet.E,
          "; ",
          "pt: ", sqrt(jet._pt2), " eta: ", rapidity(jet), " phi: ", phi(jet), ", m: ",
          m(jet), ")")
end

"""
    set_momentum!(j::PseudoJet, px, py, pz, E)

Set the momentum components and energy of a `PseudoJet` object.

# Arguments
- `j::PseudoJet`: The `PseudoJet` object to set the momentum for.
- `px`: The x-component of the momentum.
- `py`: The y-component of the momentum.
- `pz`: The z-component of the momentum.
- `E`: The energy of the particle.
"""
set_momentum!(j::PseudoJet, px, py, pz, E) = begin
    j.px = px
    j.py = py
    j.pz = pz
    j.E = E
    j._pt2 = px^2 + py^2
    j._inv_pt2 = 1.0 / j._pt2
    j._rap = _invalid_rap
    j._phi = _invalid_phi
end

"""
    _ensure_valid_rap_phi(p::PseudoJet)

Ensure that the rapidity and azimuthal angle of the PseudoJet `p` are valid. If
the azimuthal angle is invalid (used as a proxy for both variables), they are
set to a valid value using `_set_rap_phi!`.

# Arguments
- `p::PseudoJet`: The PseudoJet object to ensure valid rapidity and azimuthal
  angle for.
"""
_ensure_valid_rap_phi(p::PseudoJet) = p._phi == _invalid_phi && _set_rap_phi!(p)

"""
_set_rap_phi!(p::PseudoJet)

Set the rapidity and azimuthal angle of the PseudoJet `p`.

# Arguments
- `p::PseudoJet`: The PseudoJet object for which to set the rapidity and
  azimuthal angle.

# Description
This function calculates and sets the rapidity and azimuthal angle of the
PseudoJet `p` based on its momentum components. The rapidity is calculated in a
way that is insensitive to roundoff errors when the momentum components are
large. If the PseudoJet represents a point with infinite rapidity, a large
number is assigned to the rapidity in order to lift the degeneracy between
different zero-pt momenta.

Note - the ϕ angle is calculated in the range [0, 2π).
"""
_set_rap_phi!(p::PseudoJet) = begin
    p._phi = p._pt2 == 0.0 ? 0.0 : atan(p.py, p.px)
    if p._phi < 0.0
        p._phi += 2π
    elseif p._phi >= 2π
        p._phi -= 2π  # can happen if phi=-|eps<1e-15|?
    end

    if p.E == abs(p.pz) && iszero(p._pt2)
        # Point has infinite rapidity -- convert that into a very large
        #    number, but in such a way that different 0-pt momenta will have
        #    different rapidities (so as to lift the degeneracy between
        #                         them) [this can be relevant at parton-level]
        MaxRapHere = _MaxRap + abs(p.pz)
        p._rap = p.pz >= 0.0 ? MaxRapHere : -MaxRapHere
    else
        # get the rapidity in a way that's modestly insensitive to roundoff
        # error when things pz,E are large (actually the best we can do without
        # explicit knowledge of mass)
        effective_m2 = max(0.0, m2(p)) # force non tachyonic mass
        E_plus_pz = p.E + abs(p.pz) # the safer of p+, p-
        p._rap = 0.5 * log((p._pt2 + effective_m2) / (E_plus_pz * E_plus_pz))
        if p.pz > 0
            p._rap = -p._rap
        end
    end
    nothing
end

"""
    phi(p::PseudoJet)

Compute the ϕ angle of a `PseudoJet` object `p`.

Note this function is a wrapper for `phi_02pi(p)`.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object for which to compute the azimuthal angle.

# Returns
- The azimuthal angle of `p` in the range [0, 2π).
"""
phi(p::PseudoJet) = phi_02pi(p)

"""
    phi_02pi(p::PseudoJet)

Compute the azimuthal angle of a `PseudoJet` object `p` in the range [0, 2π).

# Arguments
- `p::PseudoJet`: The `PseudoJet` object for which to compute the azimuthal angle.

# Returns
- The azimuthal angle of `p` in the range [0, 2π).
"""
phi_02pi(p::PseudoJet) = begin
    _ensure_valid_rap_phi(p)
    return p._phi
end

"""
    rapidity(p::PseudoJet)

Compute the rapidity of a `PseudoJet` object.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object for which to compute the rapidity.

# Returns
The rapidity of the `PseudoJet` object.
"""
rapidity(p::PseudoJet) = begin
    _ensure_valid_rap_phi(p)
    return p._rap
end

"""
    pt2(p::PseudoJet)

Get the squared transverse momentum of a PseudoJet.

# Arguments
- `p::PseudoJet`: The PseudoJet object.

# Returns
- The squared transverse momentum of the PseudoJet.
"""
pt2(p::PseudoJet) = p._pt2

"""
    pt(p::PseudoJet)

Compute the scalar transverse momentum (pt) of a PseudoJet.

# Arguments
- `p::PseudoJet`: The PseudoJet object for which to compute the transverse momentum.

# Returns
- The transverse momentum (pt) of the PseudoJet.
"""
pt(p::PseudoJet) = sqrt(p._pt2)

"""
    m2(p::PseudoJet)

Calculate the invariant mass squared (m^2) of a PseudoJet.

# Arguments
- `p::PseudoJet`: The PseudoJet object for which to calculate the invariant mass squared.

# Returns
- The invariant mass squared (m^2) of the PseudoJet.
"""
m2(p::PseudoJet) = (p.E + p.pz) * (p.E - p.pz) - p._pt2

"""
    mag(p::PseudoJet)

Return the magnitude of the momentum of a `PseudoJet`, `|p|`.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object for which to compute the magnitude.

# Returns
The magnitude of the `PseudoJet` object.
"""
mag(p::PseudoJet) = sqrt(muladd(p.px, p.px, p.py^2) + p.pz^2)

"""
    CosTheta(p::PseudoJet)

Compute the cosine of the angle between the momentum vector `p` and the z-axis.

# Arguments
- `p::PseudoJet`: The PseudoJet object representing the momentum vector.

# Returns
- The cosine of the angle between `p` and the z-axis.
"""
@inline function CosTheta(p::PseudoJet)
    fZ = p.pz
    ptot = mag(p)
    return ifelse(ptot == 0.0, 1.0, fZ / ptot)
end

"""
    eta(p::PseudoJet)

Compute the pseudorapidity (η) of a PseudoJet.

# Arguments
- `p::PseudoJet`: The PseudoJet object for which to compute the pseudorapidity.

# Returns
- The pseudorapidity (η) of the PseudoJet.
"""
function eta(p::PseudoJet)
    cosTheta = CosTheta(p)
    (cosTheta^2 < 1.0) && return -0.5 * log((1.0 - cosTheta) / (1.0 + cosTheta))
    fZ = p.pz
    iszero(fZ) && return 0.0
    # Warning("PseudoRapidity","transverse momentum = 0! return +/- 10e10");
    fZ > 0.0 && return 10e10
    return -10e10
end

"""
    const η = eta

Alias for the pseudorapidity function, `eta`.
"""
const η = eta

"""
    m(p::PseudoJet)

Compute the invariant mass of a `PseudoJet` object. By convention if m^2 < 0,
then -sqrt{(-m^2)} is returned.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object for which to compute the invariant
  mass.

# Returns
The invariant mass of the `PseudoJet` object.
"""
m(p::PseudoJet) = begin
    x = m2(p)
    x < 0.0 ? -sqrt(-x) : sqrt(x)
end

"""
    mass(p::PseudoJet)

Compute the invariant mass (alias for `m(p)`).

# Arguments
- `p::PseudoJet`: The PseudoJet object for which to compute the mass.

# Returns
- The mass of the PseudoJet.
"""
mass(p::PseudoJet) = m(p)

"""Alias for `m2` function"""
const mass2 = m2

"""
    px(p::PseudoJet)

Return the x-component of the momentum of a `PseudoJet`.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object.

# Returns
- The x-component of the momentum of the `PseudoJet`.
"""
px(p::PseudoJet) = p.px

"""
    py(p::PseudoJet)

Return the y-component of the momentum of a `PseudoJet`.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object.

# Returns
- The y-component of the momentum of the `PseudoJet`.
"""
py(p::PseudoJet) = p.py

"""
    pz(p::PseudoJet)

Return the z-component of the momentum of a `PseudoJet`.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object.

# Returns
- The z-component of the momentum of the `PseudoJet`.
"""
pz(p::PseudoJet) = p.pz

"""
    energy(p::PseudoJet)

Return the energy of a `PseudoJet`.

# Arguments
- `p::PseudoJet`: The `PseudoJet` object.

# Returns
- The energy of the `PseudoJet`.
"""
energy(p::PseudoJet) = p.E

import Base.+;

"""
    +(j1::PseudoJet, j2::PseudoJet)

Addition operator for `PseudoJet` objects.

# Arguments
- `j1::PseudoJet`: The first `PseudoJet` object.
- `j2::PseudoJet`: The second `PseudoJet` object.

# Returns
A new `PseudoJet` object with the sum of the momenta and energy of `j1` and `j2`.
"""
+(j1::PseudoJet, j2::PseudoJet) = begin
    PseudoJet(j1.px + j2.px, j1.py + j2.py,
              j1.pz + j2.pz, j1.E + j2.E)
end
