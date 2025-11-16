import Accessors

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

Constructs an `EEJet` object from the given momentum components, energy, and
cluster history index.

# Details

If the default value of `cluster_hist_index=0` is used, the `EEJet` cannot be
used in a reconstruction sequence.
"""
function EEJet(px::Real, py::Real, pz::Real, E::Real; cluster_hist_index::Int = 0)
    @muladd p2 = px * px + py * py + pz * pz
    inv_p = @fastmath 1.0 / sqrt(p2)
    EEJet(px, py, pz, E, p2, inv_p, cluster_hist_index)
end

"""
    EEJet(jet::EEJet; cluster_hist_index::Int = 0)

Construct a PseudoJet from another `EEJet` object and assign given cluster index to it.
"""
function EEJet(jet::EEJet; cluster_hist_index::Int = 0)
    Accessors.@set jet._cluster_hist_index = cluster_hist_index
end

"""
    EEJet(jet::LorentzVector; cluster_hist_index::Int = 0)

Construct a EEJet from a `LorentzVector` object with optional cluster index.
"""
function EEJet(jet::LorentzVector; cluster_hist_index::Int = 0)
    EEJet(jet.x, jet.y, jet.z, jet.t; cluster_hist_index = cluster_hist_index)
end

"""
    EEJet(jet::Any; cluster_hist_index::Int = 0)

Construct a EEJet from a generic object `jet` with the given cluster index.
This functions also for `LorentzVectorCyl` objects.

# Details

This function is used to convert a generic object `jet` into an `EEJet`, for
this to work the object must have the methods `px`, `py`, `pz`, and `energy`
defined, which are used to extract the four-momentum components of the object.

The `cluster_hist_index` is optional, but needed if the `jet` is part of a
reconstruction sequence. If not provided, it defaults to `0` as an "invalid"
value.
"""
function EEJet(jet::Any; cluster_hist_index::Int = 0)
    EEJet(px(jet), py(jet), pz(jet), energy(jet);
          cluster_hist_index = cluster_hist_index)
end

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
