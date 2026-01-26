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
    EEJet(;pt::Real, rap::Real, phi::Real, m::Real = 0, cluster_hist_index::Int = 0)

Construct an EEJet from `(pt, y, ϕ, m)` with the cluster index
`cluster_hist_index`. This is not recommended, as the performance is quite
poor, but is included for completeness and to allow support for PtScheme and
Pt2Scheme.

# Details

If the (default) value of `cluster_hist_index=0` is used, the PseudoJet cannot be
used in a reconstruction sequence.
"""
function EEJet(; pt::Real, rap::Real, phi::Real, m::Real = 0,
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

    EEJet(px, py, pz, E; cluster_hist_index)
end

"""
    EEJet(jet::LorentzVector; cluster_hist_index::Int = 0)

Construct a EEJet from a `LorentzVector` object with optional cluster index.

The `cluster_hist_index` is optional, but needed if the `jet` is part of a
reconstruction sequence. If not provided, it defaults to `0` as an "invalid"
value.
"""
function EEJet(jet::LorentzVector; cluster_hist_index::Int = 0)
    EEJet(jet.x, jet.y, jet.z, jet.t; cluster_hist_index)
end

"""
    EEJet(jet::Any; cluster_hist_index::Int = 0)

Construct an EEJet from a generic object `jet` with the given cluster index.
These generic jets must implement the LorentzVectorBase interface from
`LorentzVectorBase.jl`.

# Details

This function is used to convert a generic object `jet` into a `EEJet`.
These generic jets should implement the LorentzVectorBase interface: the
`LorentzVectorBase.{px,py,pz,energy}()` methods will be called to retrieve
the four momentum components.

The `cluster_hist_index` is optional, but needed if the `jet` is part of a
reconstruction sequence. If not provided, it defaults to `0` as an "invalid"
value.
"""
function EEJet(jet::Any; cluster_hist_index::Int = 0)
    # Check that the interface is implemented
    if hasmethod(LorentzVectorBase.coordinate_system, (typeof(jet),))
        return EEJet(LorentzVectorBase.px(jet), LorentzVectorBase.py(jet),
                     LorentzVectorBase.pz(jet), LorentzVectorBase.energy(jet);
                     cluster_hist_index)
    else
        throw(ArgumentError("EEJet cannot be constructed from object of type '$(typeof(jet))'"))
    end
end

"""
    p2(eej::EEJet)

Return the squared momentum of the `EEJet` object `eej`. This accessor uses the
pre-calculated value that the struct has.
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
