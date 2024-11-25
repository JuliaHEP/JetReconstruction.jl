"""
    struct EEjet{T <: Real} <: FourMomentum

The `EEjet` struct is a 4-momentum object used for the e+e jet reconstruction
routines. Internal fields are used to track the reconstruction and to cache
values needed during the execution of the algorithm.

# Fields
- `px::T`: The x-component of the jet momentum.
- `py::T`: The y-component of the jet momentum.
- `pz::T`: The z-component of the jet momentum.
- `E::T`: The energy of the jet.
- `_cluster_hist_index::Int`: The index of the cluster histogram.
- `_p2::T`: The squared momentum of the jet.
- `_inv_p::T`: The inverse momentum of the jet.

# Type Parameters
- `T <: Real`: The type of the numerical values.
"""
mutable struct EEjet{T <: Real} <: FourMomentum
    px::T
    py::T
    pz::T
    E::T
    _p2::T
    _inv_p::T
    _cluster_hist_index::Int
end

function EEjet(px::T, py::T, pz::T, E::T, _cluster_hist_index::Integer) where {T <: Real}
    @muladd p2 = px * px + py * py + pz * pz
    inv_p = @fastmath 1.0 / sqrt(p2)
    EEjet{T}(px, py, pz, E, p2, inv_p, _cluster_hist_index)
end

EEjet(px::T, py::T, pz::T, E::T) where {T <: Real} = EEjet(px, py, pz, E, 0)

EEjet(pj::PseudoJet) = EEjet(px(pj), py(pj), pz(pj), energy(pj), cluster_hist_index(pj))

p2(eej::EEjet) = eej._p2
pt2(eej::EEjet) = eej.px^2 + eej.py^2
const kt2 = pt2
pt(eej::EEjet) = sqrt(pt2(eej))
energy(eej::EEjet) = eej.E
px(eej::EEjet) = eej.px
py(eej::EEjet) = eej.py
pz(eej::EEjet) = eej.pz
nx(eej::EEjet) = eej.px * eej._inv_p
ny(eej::EEjet) = eej.py * eej._inv_p
nz(eej::EEjet) = eej.pz * eej._inv_p
cluster_hist_index(eej::EEjet) = eej._cluster_hist_index

phi(eej::EEjet) = begin
    phi = pt2(eej) == 0.0 ? 0.0 : atan(eej.py, eej.px)
    if phi < 0.0
        phi += 2π
    elseif phi >= 2π
        phi -= 2π  # can happen if phi=-|eps<1e-15|?
    end
    phi
end

m2(eej::EEjet) = energy(eej)^2 - p2(eej)
mass(eej::EEjet) = m2(eej) < 0.0 ? -sqrt(-m2(eej)) : sqrt(m2(eej))

function rapidity(eej::EEjet)
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
function +(jet1::EEjet, jet2::EEjet)
    EEjet(jet1.px + jet2.px, jet1.py + jet2.py, jet1.pz + jet2.pz, jet1.E + jet2.E)
end

import Base.show
function show(io::IO, eej::EEjet)
    print(io, "EEjet(px: ", eej.px, " py: ", eej.py, " pz: ", eej.pz, " E: ", eej.E,
          " cluster_hist_index: ", eej._cluster_hist_index, ")")
end

"""
    mutable struct EERecoJet{T <: Real}

Optimised struct for e+e jets reconstruction, to be used with StructArrays.

# Fields
- `index::Int`: The index of the jet.
- `nni::Int`: The nearest neighbour index.
- `nndist::T`: The distance to the nearest neighbour.
- `dijdist::T`: The distance between jets.
- `nx::T`: The x-component of the jet's momentum.
- `ny::T`: The y-component of the jet's momentum.
- `nz::T`: The z-component of the jet's momentum.
- `E2p::T`: The energy raised to the power of 2p for this jet.

# Type Parameters
- `T <: Real`: The type of the numerical values.
"""
mutable struct EERecoJet{T <: Real}
    index::Int
    nni::Int
    nndist::T
    dijdist::T
    nx::T
    ny::T
    nz::T
    E2p::T
end
