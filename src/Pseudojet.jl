# Inspired by the PseudoJet class of c++ code of Fastjet (https://fastjet.fr,
#  hep-ph/0512210,  arXiv:1111.6097)
#
# Some of the implementation is taken from LorentzVectorHEP.jl, by Jerry Ling

"""
    mutable struct PseudoJet <: FourMomentum

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
- `_inv_pt2::Float64`: The inverse squared transverse momentum.
- `_rap::Float64`: The rapidity.
- `_phi::Float64`: The azimuthal angle.

"""
struct PseudoJet <: Jet
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
    _cluster_hist_index::Int
    _pt2::Float64
    _inv_pt2::Float64
    _rap::Float64
    _phi::Float64
    function PseudoJet(px, py, pz, E, cluster_hist_index)
        @muladd pt2 = px * px + py * py
        inv_pt2 = @fastmath 1.0 / pt2
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
        new(px, py, pz, E, cluster_hist_index, pt2, inv_pt2, rap, phi)
    end
end

"""
    PseudoJet(px::Real, py::Real, pz::Real, E::Real)

Constructs a PseudoJet object with the given momentum components and energy.

# Returns
A PseudoJet object with **no** cluster index.
"""
PseudoJet(px::Real, py::Real, pz::Real, E::Real) = PseudoJet(px, py, pz, E, 0)

"""
    PseudoJet(jet::PlainJet, cluster_hist_index::Int)

Constructs a PseudoJet object from a PlainJet object and a cluster history index.
"""
PseudoJet(jet::PlainJet, cluster_hist_index::Int) = PseudoJet(jet.px, jet.py, jet.pz, jet.E,
                                                              cluster_hist_index)

"""
    PseudoJet(jet::PlainJet)

Constructs a PseudoJet object from a PlainJet object with no cluster history index.
"""
PseudoJet(jet::PlainJet) = PseudoJet(jet, 0)

"""
    PseudoJet(jet::PseudoJet)

Create a copy of a PseudoJet object.
"""
PseudoJet(jet::PseudoJet) = deepcopy(jet)

"""
    PseudoJet(jet::LorentzVector, cluster_hist_index::Int)

Used to mark an invalid result in case the corresponding substructure tagging fails."""
const invalid_pseudojet = PseudoJet(0.0, 0.0, 0.0, 0.0)

import Base.isvalid
"""
    isvalid(j::PseudoJet)

Function to check whether a given `PseudoJet` object is non-zero or not.
Primarily to use for checking the validity of outputs of substructure tagging.

# Returns
- `Bool`: `true` if the `PseudoJet` object is non-zero (valid), `false` otherwise. 
"""
isvalid(j::PseudoJet) = !(j === invalid_pseudojet)

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
          "pt: ", sqrt(jet._pt2), " rapidity: ", rapidity(jet), " phi: ", phi(jet), ", m: ",
          m(jet), ")")
end

"""
    phi(p::PseudoJet)

Return the azimuthal angle, ϕ, of a `PseudoJet` object `p` in the range [0, 2π).
"""
phi(p::PseudoJet) = p._phi

"""
    rapidity(p::PseudoJet)

Compute the rapidity of a `PseudoJet` object.

# Returns
The rapidity of the `PseudoJet` object.
"""
rapidity(p::PseudoJet) = p._rap

"""
    pt2(p::PseudoJet)

Get the squared transverse momentum of a PseudoJet.

# Returns
- The squared transverse momentum of the PseudoJet.
"""
pt2(p::PseudoJet) = p._pt2

"""
    pt(p::PseudoJet)

Compute the scalar transverse momentum (pt) of a PseudoJet.

# Returns
- The transverse momentum (pt) of the PseudoJet.
"""
pt(p::PseudoJet) = sqrt(p._pt2)
