# Additional generic utility functions for jet structs
#
# Functions that create each jet type from the other need to be defined here,
# after the structs are defined.
"""
    PseudoJet(jet::EEJet; cluster_hist_index::Int=0)

Constructs a `PseudoJet` from an `EEJet`.

# Details

The `cluster_hist_index` is set to the value of the `cluster_hist_index` of the
EEJet if `0` is passed. Otherwise it is set to the value, `>0`, passed in.
"""
function PseudoJet(eej::EEJet; cluster_hist_index::Int = 0)
    cluster_hist_index = cluster_hist_index == 0 ? eej._cluster_hist_index :
                         cluster_hist_index
    PseudoJet(px(eej), py(eej), pz(eej), energy(eej);
              cluster_hist_index)
end

"""
    EEJet(jet::PseudoJet; cluster_hist_index::Int=0)

Constructs an `EEJet` from a `PseudoJet`.

# Details

The `cluster_hist_index` is set to the value of the `cluster_hist_index` of the
PseudoJet if `0` is passed. Otherwise it is set to the value, `>0`, passed in.
"""
function EEJet(jet::PseudoJet; cluster_hist_index::Int = 0)
    cluster_hist_index = cluster_hist_index == 0 ? jet._cluster_hist_index :
                         cluster_hist_index
    EEJet(px(jet), py(jet), pz(jet), energy(jet);
          cluster_hist_index)
end

# Functions to convert jets to types from other packages
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

# Utility functions for jet structs using pairs of jets
"""
    deltaR(jet1::T, jet2::T) where T <: FourMomentum

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
    deltar(jet1::T, jet2::T) where T <: FourMomentum

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
    delta_phi(jet1::T, jet2::T) where {T <: FourMomentum}

Computes the difference in azimuthal angle φ between two jets,
wrapped into the range [-π, π].

# Returns
- `δφ` as a Float64 in radians.
"""
function delta_phi(jet1::T, jet2::T) where {T <: FourMomentum}
    δϕ = phi(jet1) - phi(jet2)
    δϕ = abs(δϕ) > π ? 2π - abs(δϕ) : δϕ

    return δϕ
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
    construct_reco_jets(particles, ::Type{J}, preprocess) where {J <: FourMomentum}

Create independently owned, reconstruction-ready jets of type `J` from
`particles`.

When `preprocess` is `nothing`, inputs that already have type `J` are copied
directly and other inputs are converted to `J`. Otherwise `preprocess` is called
for every input particle.
"""
function construct_reco_jets(particles::AbstractVector{P},
                             ::Type{J},
                             preprocess) where {P, J <: FourMomentum}
    TargetNumericalType = eltype(particles[1])
    if TargetNumericalType <: Real
        TargetJetType = concretize_return_type(J, TargetNumericalType)
    else
        TargetJetType = typeof(J(particles[1]))
    end

    recombination_particles = Vector{TargetJetType}()
    sizehint!(recombination_particles, 2 * length(particles))

    return construct_reco_jets!(recombination_particles,
                                particles,
                                preprocess)
end

"""
    construct_reco_jets!(recombination_particles, particles, preprocess)

Reset and fill reusable reconstruction-jet storage. The destination must not
alias the input collection.
"""
function construct_reco_jets!(recombination_particles::Vector{J},
                              particles::AbstractVector{P},
                              preprocess) where {P, J <: FourMomentum}
    Base.mightalias(recombination_particles, particles) &&
        throw(ArgumentError("reusable jet storage must not alias the input particles"))

    N = length(particles)
    empty!(recombination_particles)
    _sizehint_for_reuse!(recombination_particles, 2 * N)

    if isnothing(preprocess)
        if P === J
            append!(recombination_particles, particles)
        else
            for (i, particle) in enumerate(particles)
                push!(recombination_particles,
                      J(particle; cluster_hist_index = i))
            end
        end
    else
        for (i, particle) in enumerate(particles)
            push!(recombination_particles,
                  preprocess(particle,
                             J;
                             cluster_hist_index = i))
        end
    end

    return recombination_particles
end
