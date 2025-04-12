"""
    abstract type FourMomentum end

Interface for composite types that includes fields px, py, py, and E
that represents the components of a four-momentum vector. All concrete
jets that are used in the package are subtypes of this type.
"""
abstract type FourMomentum end

# Define here all common functions that can be used for all jet types <: FourMomentum
import Base.+;
"""
    +(jet1::T, jet2::T) where {T <: FourMomentum}

Adds two four-momentum vectors together, returning a new jet.

# Details

This addition operation will return a jet with the cluster history index set to
0. This means that this jet cannot be used, or be part of, any clustering
history.
"""
function +(jet1::T, jet2::T) where {T <: FourMomentum}
    T(jet1.px + jet2.px, jet1.py + jet2.py, jet1.pz + jet2.pz, jet1.E + jet2.E)
end

"""
    addjets(j1::T, j2::T, cluster_hist_index::Int) where {T <: FourMomentum}

Add jets' four momenta together, returning a new jet with the specified cluster
history index.
"""
function addjets(j1::T, j2::T, cluster_hist_index::Int) where {T <: FourMomentum}
     T(j1.px+j2.px, j1.py+j2.py, j1.pz+j2.pz, j1.E+j2.E, cluster_hist_index)
end

"""
    px(j::FourMomentum)

Return the x-component of the four-momentum vector of `j`.
"""
px(j::FourMomentum) = j.px

"""
    py(j::FourMomentum)
    
Return the y-component of the four-momentum vector of `j`.
"""
py(j::FourMomentum) = j.py

"""
    pz(j::FourMomentum)
Return the z-component of the four-momentum vector of `j`.
"""
pz(j::FourMomentum) = j.pz

"""
    energy(j::FourMomentum)

Return the energy component of the four-momentum vector of `j`.
"""
energy(j::FourMomentum) = j.E

"""Alias for `energy` function"""
const E = energy

"""
    cluster_hist_index(j::FourMomentum)

Return the cluster history index of the jet `j`.
"""
cluster_hist_index(j::FourMomentum) = j._cluster_hist_index

"""
    p2(j::FourMomentum)

Return the squared momentum of the four-momentum vector of `j`.
"""
p2(j::FourMomentum) = j.px^2 + j.py^2 + j.pz^2

"""
    pt2(j::FourMomentum)

Return the squared transverse momentum of the four-momentum vector of `j`.
"""
pt2(j::FourMomentum) = j.px^2 + j.py^2

"""Alias for `pt2` function"""
const kt2 = pt2

"""
    pt(j::FourMomentum)
Return the momentum of the four-momentum vector of `j`.
"""
pt(j::FourMomentum) = sqrt(pt2(j))

"""
    mass2(j::FourMomentum)

Return the invariant mass squared of the four-momentum vector `j`.
"""
mass2(j::FourMomentum) = energy(j)^2 - p2(j)

"""Alias for `mass2` function"""
const m2 = mass2

"""
    mass(j::FourMomentum)

Return the invariant mass of a four momentum `j`. By convention if ``m^2 < 0``,
then ``-sqrt{(-m^2)}`` is returned.
"""
mass(j::FourMomentum) = m2(j) < 0.0 ? -sqrt(-m2(j)) : sqrt(m2(j))

"""Alias for `mass` function"""
const m = mass

"""
    phi(j::FourMomentum)

Return the azimuthal angle, ϕ, of the four momentum `j` in the range [0, 2π).
"""
phi(j::FourMomentum) = begin
    phi = pt2(j) == 0.0 ? 0.0 : atan(j.py, j.px)
    if phi < 0.0
        phi += 2π
    end
    phi
end

"""Used to protect against parton-level events where pt can be zero
for some partons, giving rapidity=infinity. KtJet fails in those cases."""
const _MaxRap = 1e5

"""
    rapidity(j::FourMomentum)

Return the rapidity of the four momentum `j`.
"""
function rapidity(j::FourMomentum)
    if energy(j) == abs(pz(j)) && iszero(pt2(j))
        MaxRapHere = _MaxRap + abs(pz(j))
        rapidity = (pz(j) >= 0.0) ? MaxRapHere : -MaxRapHere
    else
        # get the rapidity in a way that's modestly insensitive to roundoff
        # error when things pz,E are large (actually the best we can do without
        # explicit knowledge of mass)
        effective_m2 = max(0.0, m2(j)) # force non tachyonic mass
        E_plus_pz = energy(j) + abs(pz(j)) # the safer of p+, p-
        rapidity = 0.5 * log((pt2(j) + effective_m2) / (E_plus_pz * E_plus_pz))
        if pz(j) > 0
            rapidity = -rapidity
        end
    end
    rapidity
end
