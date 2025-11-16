"""
    abstract type FourMomentum end

Interface for composite types that includes fields px, py, py, and E
that represents the components of a four-momentum vector. All concrete
jets that are used in the package are subtypes of this type.
"""
abstract type FourMomentum end

# Define here all common functions that can be used for all jet types <: FourMomentum
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
    p(j::FourMomentum)

Return the momentum of the four-momentum vector of `j`.
"""
p(j::FourMomentum) = sqrt(p2(j))

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

# Alternative name for [0, 2π) range
const phi02pi = phi

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

"""
    mag(jet::T) where {T <: FourMomentum}

Return the magnitude of the momentum of a jet, `|p|`.

# Returns
The magnitude of the jet.
"""
mag(jet::T) where {T <: FourMomentum} = sqrt(muladd(jet.px, jet.px, jet.py^2) + jet.pz^2)

"""
    CosTheta(jet::T) where {T <: FourMomentum}

Compute the cosine of the angle between the momentum vector of `jet` and the z-axis.

# Returns
- The cosine of the angle between the jet and the z-axis.
"""
@inline function CosTheta(jet::T) where {T <: FourMomentum}
    fZ = jet.pz
    ptot = mag(jet)
    return ifelse(ptot == 0.0, 1.0, fZ / ptot)
end

"""
    eta(jet::T) where {T <: FourMomentum}

Compute the pseudorapidity (η) of a jet.

# Returns
- The pseudorapidity (η) of the jet.
"""
function eta(jet::T) where {T <: FourMomentum}
    cosTheta = CosTheta(jet)
    (cosTheta^2 < 1.0) && return -0.5 * log((1.0 - cosTheta) / (1.0 + cosTheta))
    fZ = jet.pz
    iszero(fZ) && return 0.0
    # Warning("PseudoRapidity","transverse momentum = 0! return +/- 10e10");
    fZ > 0.0 ? _MaxRap : -_MaxRap
end

"""
    const η = eta

Alias for the pseudorapidity function, `eta`.
"""
const η = eta

import Base.+;
"""
    +(jet1::T, jet2::T) where {T <: FourMomentum}

Adds two four-momentum vectors together, returning a new jet.

# Details

This addition operation will return a jet with the cluster history index set to
0. *This means that this jet cannot be used, or be part of, any clustering
history.*
"""
function +(jet1::T, jet2::T) where {T <: FourMomentum}
    T(jet1.px + jet2.px, jet1.py + jet2.py, jet1.pz + jet2.pz, jet1.E + jet2.E)
end

"""
    addjets(jet1::T, jet2::T; cluster_hist_index::Int) where {T <: FourMomentum}

Add jets' four momenta together, returning a new jet of type `T` with the
specified cluster history index.

# Details

This method is also known as the `E_scheme` in Fastjet.
"""
function addjets(jet1::T, jet2::T; cluster_hist_index::Int = 0) where {T <: FourMomentum}
    T(px(jet1) + px(jet2), py(jet1) + py(jet2), pz(jet1) + pz(jet2),
      energy(jet1) + energy(jet2), cluster_hist_index = cluster_hist_index)
end

"""
    const addjets_escheme = addjets

Alias for the `addjets` function, which is the default jet recombination scheme.
"""
const addjets_escheme = addjets

"""
    addjets_ptscheme(jet1::T, jet2::T, cluster_hist_index::Int) where {T <: FourMomentum}

Use the massless ``p_T`` scheme for combining two jets, setting the appropriate
cluster history index for the new jet.
"""
function addjets_ptscheme(jet1::T, jet2::T;
                          cluster_hist_index::Int = 0) where {T <: FourMomentum}
    scale1 = pt(jet1)
    scale2 = pt(jet2)
    _addjets_with_scale(scale1, scale2, jet1, jet2, cluster_hist_index)::T
end

"""
    preprocess_escheme(jet::T, ::Type{OutputT};
                                cluster_hist_index::Int = 0) -> OutputT where {T <: FourMomentum,
                                                                               OutputT <: FourMomentum}

Jet preprocessor for the E-scheme, simply copying the four-momentum and assigning the cluster
history index.
"""
function preprocess_escheme(jet::T, ::Type{OutputT};
                            cluster_hist_index::Int = 0) where {T,
                                                                OutputT <: FourMomentum}
    OutputT(jet; cluster_hist_index = cluster_hist_index)
end

"""
    preprocess_escheme(jet::T; cluster_hist_index::Int = 0) -> T where {T <: FourMomentum}

Jet preprocessor for the E-scheme, simply copying the four-momentum and assigning the cluster
history index.
"""
function preprocess_escheme(jet::T;
                            cluster_hist_index::Int = 0) where {T <: FourMomentum}
    preprocess_escheme(jet, T; cluster_hist_index = cluster_hist_index)
end

"""
    preprocess_ptscheme(jet::T, ::Type{OutputT};
                             cluster_hist_index::Int = 0) -> OutputT where {T <: FourMomentum,
                                                                            OutputT <: FourMomentum}

Jet preprocessor for the massless ``p_T`` schemes, resetting the energy of the
jet to be equal to the 3-momentum of the input jet.
"""
function preprocess_ptscheme(jet::T, ::Type{OutputT};
                             cluster_hist_index::Int = 0) where {T <: FourMomentum,
                                                                 OutputT <: FourMomentum}
    OutputT(px(jet), py(jet), pz(jet), p(jet); cluster_hist_index = cluster_hist_index)
end

"""
    preprocess_ptscheme(jet::T; cluster_hist_index::Int = 0) -> T where {T <: FourMomentum}

Jet preprocessor for the massless ``p_T`` schemes, resetting the energy of the
jet to be equal to the 3-momentum of the input jet.
"""
function preprocess_ptscheme(jet::T;
                             cluster_hist_index::Int = 0) where {T <: FourMomentum}
    preprocess_ptscheme(jet, T; cluster_hist_index = cluster_hist_index)
end

"""
    preprocess_ptscheme(particle::Union{LorentzVector, LorentzVectorCyl},
                             ::Type{OutputT}= PseudoJet,;
                             cluster_hist_index::Int = 0) -> OutputT where {OutputT <: FourMomentum}

Jet preprocessor for the massless ``p_T`` schemes, resetting the energy of the
jet to be equal to the 3-momentum of the input jet (generic particle type).

# Details

This function is used to convert a particle of type `LorentzVector` or
`LorentzVectorCyl` into a `jet_type` object, which is a subtype of
`FourMomentum`. (This is a work around until `LorentzVectorBase` can be used,
which will make the accessors uniform.)
"""
function preprocess_ptscheme(particle::Union{LorentzVector, LorentzVectorCyl},
                             ::Type{OutputT} = PseudoJet;
                             cluster_hist_index::Int = 0) where {OutputT <: FourMomentum}
    OutputT(px(particle), py(particle), pz(particle), mag(particle);
            cluster_hist_index = cluster_hist_index)
end

"""
    addjets_pt2scheme(jet1::T, jet2::T, cluster_hist_index::Int) where {T <: FourMomentum}

Use the massless ``p_T^2`` scheme for combining two jets, setting the
appropriate cluster history index for the new jet.
"""
function addjets_pt2scheme(jet1::T, jet2::T;
                           cluster_hist_index::Int) where {T <: FourMomentum}
    scale1 = pt2(jet1)
    scale2 = pt2(jet2)
    _addjets_with_scale(scale1, scale2, jet1, jet2, cluster_hist_index)::T
end

"""
    const preprocess_pt2scheme = preprocess_ptscheme

Preprocessing for `p_T` and `p_T^2` schemes are identical.
"""
const preprocess_pt2scheme = preprocess_ptscheme

"""
    _addjets_with_scale(scale1::Real, scale2::Real, jet1::T, jet2::T, cluster_hist_index::Int) where {T <: FourMomentum}

Combine two jets as massless objects using the given scale factors for each jet,
and return the new jet with the cluster history index set.
"""
function _addjets_with_scale(scale1::Real, scale2::Real, jet1::T, jet2::T,
                             cluster_hist_index::Int) where {T <: FourMomentum}
    new_pt = pt(jet1) + pt(jet2)
    new_rap = (scale1 * rapidity(jet1) + scale2 * rapidity(jet2)) / (scale1 + scale2)
    phi_wrap = 0.0
    if phi(jet1) - phi(jet2) > π
        phi_wrap = 2π
    elseif phi(jet1) - phi(jet2) < -π
        phi_wrap = -2π
    end
    @assert -π < phi(jet1) - (phi(jet2) + phi_wrap) < π
    new_phi = (scale1 * phi(jet1) + scale2 * (phi(jet2) + phi_wrap)) / (scale1 + scale2)

    # Now create new jet from pt, y and phi... implicitly assume that mass=0!
    T(pt = new_pt, rap = new_rap, phi = new_phi; cluster_hist_index = cluster_hist_index)
end

import Base.show
"""
    show(io::IO, jet::FourMomentum)

Print core information about the four-momentum vector of `jet` to the given IO
stream.
"""
function show(io::IO, jet::FourMomentum)
    print(io, "$(typeof(jet))(px: ", px(jet), " py: ", py(jet), " pz: ", pz(jet), " E: ",
          E(jet),
          " cluster_hist_index: ", cluster_hist_index(jet), ")")
end
