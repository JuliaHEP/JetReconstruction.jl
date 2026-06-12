# Additional generic utility functions for jet structs
#
# Functions that create each jet type from the other need to be defined here,
# after the structs are defined.
"""
    PseudoJet(eej::EEJet; cluster_hist_index::Int = 0)

Constructs a `PseudoJet` from an `EEJet`. The underlying type `T` is preserved.

# Details

The `cluster_hist_index` is set to the value of the `cluster_hist_index` of the
`EEJet` if `0` is passed. Otherwise it is set to the value, `>0`, passed in.
"""
function PseudoJet(eej::EEJet; cluster_hist_index::Int = 0)
    cluster_hist_index = cluster_hist_index == 0 ? eej._cluster_hist_index :
                         cluster_hist_index
    PseudoJet(px(eej), py(eej), pz(eej), energy(eej);
              cluster_hist_index)
end

"""
    EEJet(jet::PseudoJet; cluster_hist_index::Int = 0)

Constructs an `EEJet` from a `PseudoJet`. The underlying type `T` is preserved.

# Details

The `cluster_hist_index` is set to the value of the `cluster_hist_index` of the
`PseudoJet` if `0` is passed. Otherwise it is set to the value, `>0`, passed in.
"""
function EEJet(jet::PseudoJet; cluster_hist_index::Int = 0)
    cluster_hist_index = cluster_hist_index == 0 ? jet._cluster_hist_index :
                         cluster_hist_index
    EEJet(px(jet), py(jet), pz(jet), energy(jet);
          cluster_hist_index)
end

# Functions to convert jets to types from other packages
"""
    lorentzvector_cyl(jet::FourMomentum)

Return a cylindrical `LorentzVectorCyl` from a jet.
"""
function lorentzvector_cyl(jet::T) where {T <: FourMomentum}
    return LorentzVectorHEP.fromPtEtaPhiE(JetReconstruction.pt(jet),
                                          JetReconstruction.eta(jet),
                                          JetReconstruction.phi(jet),
                                          JetReconstruction.energy(jet))
end

"""
    lorentzvector(jet::FourMomentum)

Return a cartesian `LorentzVector` from a jet.
"""
function lorentzvector(jet::T) where {T <: FourMomentum}
    return LorentzVector(JetReconstruction.energy(jet),
                         JetReconstruction.px(jet),
                         JetReconstruction.py(jet),
                         JetReconstruction.pz(jet))
end

# Utility functions for jet structs using pairs of jets
"""
    deltaR(jet1::FourMomentum, jet2::FourMomentum)

Function to calculate the distance in the y-ϕ plane between two jets `jet1` and
`jet2` (that is using *rapidity* and *azimuthal angle*).

# Returns
- The Euclidean distance in the y-ϕ plane for the two jets.
"""
function deltaR(jet1::T, jet2::T) where {T <: FourMomentum}
    δy = rapidity(jet1) - rapidity(jet2)
    δϕ = delta_phi(jet1, jet2)

    return sqrt(δy^2 + δϕ^2)
end

"""
    deltar(jet1::FourMomentum, jet2::FourMomentum)

Function to calculate the distance in the η-ϕ plane between two jets `jet1` and
`jet2` (that is, using the *pseudorapidity* and *azimuthal angle*).

# Returns
- The Euclidean distance in the η-ϕ plane for the two jets.
"""
function deltar(jet1::T, jet2::T) where {T <: FourMomentum}
    δη = eta(jet1) - eta(jet2)
    δϕ = delta_phi(jet1, jet2)

    return sqrt(δη^2 + δϕ^2)
end

"""
    delta_phi(jet1::FourMomentum, jet2::FourMomentum)

Computes the difference in azimuthal angle φ between two jets,
wrapped into the range [-π, π].

"""
function delta_phi(jet1::T, jet2::T) where {T <: FourMomentum}
    δϕ = phi(jet1) - phi(jet2)
    δϕ = abs(δϕ) > π ? 2π - abs(δϕ) : δϕ

    return δϕ
end

"""
    pt_fraction(jet1::FourMomentum, jet2::FourMomentum)

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
    kt_scale(jet1::FourMomentum, jet2::FourMomentum)

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
    construct_reco_jets(particles, ::Type{J}, preprocess) where {J <: FourMomentum}

Generate reconstruction ready jets from the set of input particles.

# Details

The initial "jets" ready for reconstruction will be created from `particles`. The
type of the jets is given by `J` so that different `FourMomentum` types can be
used. The function `preprocess` will be used to massage the input particles, unless
`nothing` is passed, in which case this is skipped.
"""
function construct_reco_jets(particles::P, ::Type{J},
                             preprocess) where {P, J <: FourMomentum}
    # Decide what our return type will be, based on what we find
    # in `particles`
    # See if `particles` sensibly supports `eltype()`, otherwise
    # apply a hack with an explicit conversion
    TargetNumericalType = eltype(particles[1])
    if TargetNumericalType <: Real
        TargetJetType = J{TargetNumericalType}
    else
        TargetJetType = typeof(J(particles[1]))
    end

    if isnothing(preprocess)
        if P == J
            # If we don't have a preprocessor, we just need to copy the inputs
            # when the type matches the target
            recombination_particles = copy(particles)
            sizehint!(recombination_particles, length(particles) * 2)
        else
            # We assume a constructor for J that can ingest the appropriate
            # type of particle
            recombination_particles = Vector{TargetJetType}(undef, 0)
            sizehint!(recombination_particles, length(particles) * 2)
            for (i, particle) in enumerate(particles)
                push!(recombination_particles, J(particle; cluster_hist_index = i))
            end
        end
    else
        # We have a preprocessor function that we need to call to modify the
        # input particles
        recombination_particles = Vector{TargetJetType}(undef, 0)
        sizehint!(recombination_particles, length(particles) * 2)
        for (i, particle) in enumerate(particles)
            push!(recombination_particles,
                  preprocess(particle, TargetJetType, cluster_hist_index = i))
        end
    end
    recombination_particles
end
