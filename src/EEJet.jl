"""
    struct EEJet <: FourMomentum

The `EEJet` struct is a 4-momentum object used for the e+e jet reconstruction routines.

# Fields
- `px::Float64`: The x-component of the jet momentum.
- `py::Float64`: The y-component of the jet momentum.
- `pz::Float64`: The z-component of the jet momentum.
- `E::Float64`: The energy of the jet.
- `_cluster_hist_index::Int`: The index of the cluster histogram.
- `_p2::Float64`: The squared momentum of the jet.
- `_inv_p::Float64`: The inverse momentum of the jet.
"""
struct EEJet <: FourMomentum
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
    _p2::Float64
    _inv_p::Float64
    _cluster_hist_index::Int
end

"""
    EEJet(px::Real, py::Real, pz::Real, E::Real, cluster_hist_index::Int)

Constructs an `EEJet` object from the given momentum components, energy, and cluster history index.
"""
function EEJet(px::Real, py::Real, pz::Real, E::Real, cluster_hist_index::Int)
    @muladd p2 = px * px + py * py + pz * pz
    inv_p = @fastmath 1.0 / sqrt(p2)
    EEJet(px, py, pz, E, p2, inv_p, cluster_hist_index)
end

"""
    EEJet(px::Real, py::Real, pz::Real, E::Real)

Constructs an `EEJet` object from the given momentum components and energy,
but with no cluster history index (which is set to 0).
"""
EEJet(px::Real, py::Real, pz::Real, E::Real) = EEJet(px, py, pz, E, 0)

"""
    p2(eej::EEJet)

Return the squared momentum of the `EEJet` object `eej`.
"""
p2(eej::EEJet) = eej._p2

"""
    nx(eej::EEJet)

Return the x-component of the unit vector aligned with the momentum of `eej`
"""
nx(eej::EEJet) = eej.px * eej._inv_p

"""
    ny(eej::EEJet)

Return the y-component of the unit vector aligned with the momentum of `eej`
"""
ny(eej::EEJet) = eej.py * eej._inv_p

"""
    nz(eej::EEJet)

Return the z-component of the unit vector aligned with the momentum of `eej`
"""
nz(eej::EEJet) = eej.pz * eej._inv_p

# Optimised reconstruction struct for e+e jets

mutable struct EERecoJet
    index::Int
    nni::Int
    nndist::Float64
    dijdist::Float64
    nx::Float64
    ny::Float64
    nz::Float64
    E2p::Float64
end
