# Substructure specific functions for jet grooming and filtering

"""
    recluster(jet, clusterseq; R = 1.0, algorithm = JetAlgorithm.CA) -> ClusterSequence

Reclusters the constituents of a given jet `jet` with a different clustering algorithm `algorithm` and different jet radius `R`.

# Arguments
- `jet::PseudoJet`: The jet whose constituents are to be reclustered.
- `clusterseq::ClusterSequence{PseudoJet}`: The cluster sequence from which the original jet is obtained.
- `R = 1.0`: The new jet radius.
- `algorithm::JetAlgorithm.Algorithm = JetAlgorithm.CA`: The new clustering method.

# Returns
- `ClusterSequence`: The new cluster sequence.
"""
function recluster(jet::PseudoJet, clusterseq::ClusterSequence{PseudoJet}; R = 1.0,
                   algorithm::JetAlgorithm.Algorithm = JetAlgorithm.CA)
    cons = constituents(jet, clusterseq)
    new_clusterseq = jet_reconstruct(cons; p = nothing, R = R, algorithm = algorithm,
                                     strategy = RecoStrategy.Best)

    return new_clusterseq
end

"""
    mass_drop(jet, clusterseq; mu, y) -> PseudoJet

Identifies subjets in a jet that pass the mass drop tagging condition.
The method stops at the first jet satisfying the mass and distance thresholds.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance representing the jet to tag.
- `clusterseq::ClusterSequence{PseudoJet}`: ClusterSequence with jet clustering history.
- `mu::Real`: Maximum allowed mass ratio for a jet to pass tagging.
- `y::Real`: Minimum kT distance threshold for parent separation.

# Returns:
- `PseudoJet`: The jet (or subjet) satisfying the mass drop conditions, if tagging is successful, otherwise `invalid_pseudojet` object
"""
function mass_drop(jet::PseudoJet, clusterseq::ClusterSequence{PseudoJet}; mu::Real,
                   y::Real)
    all_jets = clusterseq.jets
    hist = clusterseq.history

    while true
        parent1, parent2 = parent_jets(jet, clusterseq)

        if !isnothing(parent1) && !(isnothing(parent2))
            if m2(parent1) < m2(parent2)
                parent1, parent2 = parent2, parent1
            end

            pt1 = pt(parent1)
            pt2 = pt(parent2)

            if m2(parent1) < m2(jet) * mu^2 &&
               (min(pt1, pt2) * deltaR(parent1, parent2))^2 > y * m2(jet)
                return jet
            else
                jet = parent1
            end

        else
            return invalid_pseudojet
        end
    end
end

"""
    soft_drop(jet, clusterseq; zcut, beta, radius = 1.0) -> PseudoJet

Applies soft-drop grooming to remove soft, wide-angle radiation from jets.
This function reclusters the jet and iteratively checks the soft-drop condition on subjets.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance to groom.
- `clusterseq::ClusterSequence{PseudoJet}`: ClusterSequence containing jet history.
- `zcut::Real`: Minimum allowed energy fraction for subjets.
- `beta::Real`: Angular exponent controlling soft radiation suppression.
- `radius::Real`: The new radius that will be used to recluster the
  components of the jet, by default set to 1.0.

# Returns:
- `PseudoJet`: Groomed jet or `invalid_pseudojet` object if grooming fails.
"""
function soft_drop(jet::PseudoJet, clusterseq::ClusterSequence{PseudoJet}; zcut::Real,
                   beta::Real, radius::Real = 1.0)
    new_clusterseq = recluster(jet, clusterseq; R = radius, algorithm = JetAlgorithm.CA)
    new_jet = sort!(inclusive_jets(new_clusterseq; T = PseudoJet), by = pt2, rev = true)[1]

    all_jets = new_clusterseq.jets
    hist = new_clusterseq.history

    while true
        parent1, parent2 = parent_jets(new_jet, new_clusterseq)

        if !isnothing(parent1) && !(isnothing(parent2))
            if m2(parent1) < m2(parent2)
                parent1, parent2 = parent2, parent1
            end

            pt1 = pt(parent1)
            pt2 = pt(parent2)

            if min(pt1, pt2) / (pt1 + pt2) >
               zcut * (deltaR(parent1, parent2) / radius)^beta
                return new_jet
            else
                new_jet = parent1
            end

        else
            return invalid_pseudojet
        end
    end
end

"""
    jet_filtering(jet, clusterseq; radius, hardest_jets) -> PseudoJet

Filters a jet to retain only the hardest subjets based on a specified radius and number.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance representing the jet to filter.
- `clusterseq::ClusterSequence{PseudoJet}`: ClusterSequence containing jet history.
- `radius::Real`: Radius parameter to recluster subjets.
- `hardest_jets::Integer`: Number of hardest jets to retain in the filtered result.

# Returns:
- `PseudoJet`: Filtered jet composed of the hardest subjets.
"""
function jet_filtering(jet::PseudoJet, clusterseq::ClusterSequence{PseudoJet}; radius::Real,
                       hardest_jets::Integer)
    new_clusterseq = recluster(jet, clusterseq; R = radius, algorithm = JetAlgorithm.CA)
    reclustered = sort!(inclusive_jets(new_clusterseq; T = PseudoJet), by = pt2, rev = true)

    n = length(reclustered) <= hardest_jets ? length(reclustered) : hardest_jets
    hard = reclustered[1:n]

    filtered = length(hard) != 0 ? foldl(+, hard) : invalid_pseudojet
    return filtered
end

"""
    jet_trimming(jet, clusterseq; radius, fraction, recluster_method) -> PseudoJet

Trims a jet by removing subjets with transverse momentum below a specified fraction.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance representing the jet to trim.
- `clusterseq::ClusterSequence{PseudoJet}`: ClusterSequence containing jet history.
- `radius::Real`: Radius used for reclustering in trimming.
- `fraction::Real`: Minimum momentum fraction for retained subjets.
- `recluster_method::JetAlgorithm.Algorithm`: Method identifier for reclustering.

# Returns:
- `PseudoJet`: Trimmed jet composed of retained subjets.
"""
function jet_trimming(jet::PseudoJet, clusterseq::ClusterSequence{PseudoJet}; radius::Real,
                      fraction::Real, recluster_method::JetAlgorithm.Algorithm)
    frac2 = fraction^2

    new_clusterseq = recluster(jet, clusterseq; R = radius, algorithm = recluster_method)
    reclustered = sort!(inclusive_jets(new_clusterseq; T = PseudoJet), by = pt2, rev = true)

    hard = Vector{PseudoJet}(undef, 0)
    for item in reclustered
        if pt2(item) >= frac2 * pt2(jet)
            push!(hard, item)
        end
    end

    trimmed = length(hard) != 0 ? foldl(+, hard) : invalid_pseudojet
    return trimmed
end
