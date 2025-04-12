# Generic utility functions for jet structs

# Functions that create each jet type from the other need to be defined here,
# after the structs are defined.
"""
    PseudoJet(jet::EEJet)

Constructs a `PseudoJet` object from an `EEJet` object, with the same four
momentum and cluster history index.
"""
function PseudoJet(eej::EEJet)
    PseudoJet(px(eej), py(eej), pz(eej), energy(eej), cluster_hist_index(eej))
end

"""
    EEJet(jet::PseudoJet)

Constructs an `EEJet` object from a `PseudoJet` object using the same cluster
history index from the `PseudoJet`.
"""
EEJet(jet::PseudoJet) = EEJet(px(jet), py(jet), pz(jet), energy(jet),
                              cluster_hist_index(jet))

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
    fZ > 0.0 && return 10e10
    return -10e10
end

"""
    const η = eta

Alias for the pseudorapidity function, `eta`.
"""
const η = eta

"""
    deltaR(jet1::T, jet2::T) where T <: FourMomentum

Function to calculate the distance in the y-ϕ plane between two jets `jet1` and
`jet2` (that is using *rapidity* and *azimuthal angle*).

# Returns
- The Euclidean distance in the y-ϕ plane for the two jets.
"""
function deltaR(jet1::T, jet2::T) where {T <: FourMomentum}
    δy = rapidity(jet1) - rapidity(jet2)
    δϕ = phi(jet1) - phi(jet2)
    δϕ = abs(δϕ) > π ? 2π - abs(δϕ) : δϕ

    return sqrt(δy^2 + δϕ^2)
end

"""
    deltar(jet1::T, jet2::T) where T <: FourMomentum

Function to calculate the distance in the η-ϕ plane between two jets `jet1` and
`jet2` (that is, using the *pseudorapidity* and *azimuthal angle*).

# Returns
- The Euclidean distance in the η-ϕ plane for the two jets.
"""
function deltar(jet1::T, jet2::T) where {T <: FourMomentum}
    δη = eta(jet1) - eta(jet2)
    δϕ = phi(jet1) - phi(jet2)
    δϕ = abs(δϕ) > π ? 2π - abs(δϕ) : δϕ

    return sqrt(δη^2 + δϕ^2)
end

"""
    pt_fraction(jet1::T, jet2::T) where T <: FourMomentum

Computes the transverse momentum fraction of the softer of two jets.

# Returns
- The transverse momentum fraction of the softer of the two jets.
"""
function pt_fraction(jet1::T, jet2::T) where {T <: FourMomentum}
    pt1 = JetReconstruction.pt(jet1)
    pt2 = JetReconstruction.pt(jet2)
    return min(pt1, pt2) / (pt1 + pt2)
end

"""
    kt_scale(jet1::T, jet2::T) where {T <: FourMomentum}

Computes the transverse momentum scale as the product of the minimum pt and 
the angular separation in the η-ϕ plane (using *pseudorapidity*).

# Returns
- The transverse momentum scale of the two jets.
"""
function kt_scale(jet1::T, jet2::T) where {T <: FourMomentum}
    pt1 = JetReconstruction.pt(jet1)
    pt2 = JetReconstruction.pt(jet2)
    return min(pt1, pt2) * deltar(jet1, jet2)
end

"""
    lorentzvector_cyl(jet::T) where T <: FourMomentum -> LorentzVectorHEP.LorentzVectorCyl

Return a cylindrical `LorentzVectorCyl` from a jet.
"""
function lorentzvector_cyl(jet::T) where {T <: FourMomentum}
    return LorentzVectorHEP.fromPtEtaPhiE(JetReconstruction.pt(jet),
                                          JetReconstruction.eta(jet),
                                          JetReconstruction.phi(jet),
                                          JetReconstruction.energy(jet))
end

"""
    lorentzvector(jet::T) where {T <: FourMomentum} ->  -> LorentzVector

Return a cartesian `LorentzVector` from a jet.
"""
function lorentzvector(jet::T) where {T <: FourMomentum}
    return LorentzVector(JetReconstruction.energy(jet),
                         JetReconstruction.px(jet),
                         JetReconstruction.py(jet),
                         JetReconstruction.pz(jet))
end
