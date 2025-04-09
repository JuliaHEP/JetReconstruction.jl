"""
    mutable struct EEJet

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

function EEJet(px::Real, py::Real, pz::Real, E::Real, cluster_hist_index::Int)
    @muladd p2 = px * px + py * py + pz * pz
    inv_p = @fastmath 1.0 / sqrt(p2)
    EEJet(px, py, pz, E, p2, inv_p, cluster_hist_index)
end

EEJet(px::Real, py::Real, pz::Real, E::Real) = EEJet(px, py, pz, E, 0)

EEJet(pj::PseudoJet) = EEJet(px(pj), py(pj), pz(pj), energy(pj), cluster_hist_index(pj))

"""
Simple jet structure without cluster history
"""
struct PlainJet <: FourMomentum
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
end

EEJet(jet::PlainJet, cluster_hist_index) = EEJet(jet.px, jet.py, jet.pz, jet.E, cluster_hist_index)
EEJet(jet::LorentzVector, cluster_hist_index) = EEJet(px(jet), py(jet), pz(jet), energy(jet), cluster_hist_index)
EEJet(jet::LorentzVectorCyl, cluster_hist_index) = EEJet(px(jet), py(jet), pz(jet), energy(jet), cluster_hist_index)

p2(eej::EEJet) = eej._p2
pt2(eej::EEJet) = eej.px^2 + eej.py^2
const kt2 = pt2
pt(eej::EEJet) = sqrt(pt2(eej))
energy(eej::EEJet) = eej.E
px(eej::EEJet) = eej.px
py(eej::EEJet) = eej.py
pz(eej::EEJet) = eej.pz
nx(eej::EEJet) = eej.px * eej._inv_p
ny(eej::EEJet) = eej.py * eej._inv_p
nz(eej::EEJet) = eej.pz * eej._inv_p
cluster_hist_index(eej::EEJet) = eej._cluster_hist_index

phi(eej::EEJet) = begin
    phi = pt2(eej) == 0.0 ? 0.0 : atan(eej.py, eej.px)
    if phi < 0.0
        phi += 2π
    end
    phi
end

m2(eej::EEJet) = energy(eej)^2 - p2(eej)
mass(eej::EEJet) = m2(eej) < 0.0 ? -sqrt(-m2(eej)) : sqrt(m2(eej))

function rapidity(eej::EEJet)
    if energy(eej) == abs(pz(eej)) && iszero(pt2(eej))
        MaxRapHere = _MaxRap + abs(pz(eej))
        rap = (pz(eej) >= 0.0) ? MaxRapHere : -MaxRapHere
    else
        # get the rapidity in a way that's modestly insensitive to roundoff
        # error when things pz,E are large (actually the best we can do without
        # explicit knowledge of mass)
        effective_m2 = max(0.0, m2(eej)) # force non tachyonic mass
        E_plus_pz = energy(eej) + abs(pz(eej)) # the safer of p+, p-
        rapidity = 0.5 * log((pt2(eej) + effective_m2) / (E_plus_pz * E_plus_pz))
        if pz(eej) > 0
            rapidity = -rapidity
        end
    end
    rapidity
end

import Base.+;
function +(jet1::EEJet, jet2::EEJet)
    PlainJet(jet1.px + jet2.px, jet1.py + jet2.py, jet1.pz + jet2.pz, jet1.E + jet2.E)
    # LorentzVector{Float64}(jet1.E + jet2.E, jet1.px + jet2.px, jet1.py + jet2.py, jet1.pz + jet2.pz)
end

import Base.show
function show(io::IO, eej::EEJet)
    print(io, "EEJet(px: ", eej.px, " py: ", eej.py, " pz: ", eej.pz, " E: ", eej.E,
          " cluster_hist_index: ", eej._cluster_hist_index, ")")
end

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
