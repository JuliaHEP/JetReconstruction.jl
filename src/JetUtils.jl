# Generic utility functions for jet structs

"""
    deltaR(jet1::T, jet2::T) where T <: FourMomentum

Function to calculate the distance in the y-ϕ plane between two jets `jet1` and
`jet2` (that is using *rapidity* and *azimuthal angle*).

# Arguments
- `jet1::T`: The first jet.
- `jet2::T`: The second jet.

# Returns
- The Euclidean distance in the y-ϕ plane for the two jets.
"""
function deltaR(jet1::T, jet2::T) where T <: FourMomentum
    y1, ϕ1 = rapidity(jet1), phi(jet1)
    y2, ϕ2 = rapidity(jet2), phi(jet2)

    δy = y1 - y2
    δϕ = ϕ1 - ϕ2
    δϕ = abs(δϕ) > π ? 2π - abs(δϕ) : δϕ

    return sqrt(δy^2 + δϕ^2)
end

"""
    deltar(jet1::T, jet2::T) where T <: FourMomentum

Function to calculate the distance in the η-ϕ plane between two jets `jet1` and
`jet2` (that is, using the *pseudorapidity* and *azimuthal angle*).

# Arguments
- `jet1::T`: The first jet.
- `jet2::T`: The second jet.

# Returns
- The Euclidean distance in the η-ϕ plane for the two jets.
"""
function deltar(jet1::T, jet2::T) where T <: FourMomentum
    η1, ϕ1 = eta(jet1), phi(jet1)
    η2, ϕ2 = eta(jet2), phi(jet2)
ϕ
    δy = η1 - η2
    δϕ = ϕ1 - ϕ2
    δϕ = abs(δϕ) > π ? 2π - abs(δϕ) : δϕ

    return sqrt(δy^2 + δϕ^2)
end

"""
    momentum_fraction(jet1::T, jet2::T) where T <: FourMomentum

Computes the momentum fraction of the softer of two jets
"""
function momentum_fraction(jet1::T, jet2::T) where T <: FourMomentum
    pt1 = JetReconstruction.pt(jet1)
    pt2 = JetReconstruction.pt(jet2)
    return min(pt1, pt2) / (pt1 + pt2)
end

"""
    kt_scale(pj1::T, pj2::T) where T <: FourMomentum

Computes the transverse momentum scale as the product of the minimum pt and 
the angular separation computed via LorentzVectorHEP's deltar.
"""
function kt_scale(jet1::T, jet2::T) where T <: FourMomentum
    pt1 = JetReconstruction.pt(jet1)
    pt2 = JetReconstruction.pt(jet2)
    return min(pt1, pt2) * deltar(jet1, jet2)
end
