"""
    deltaR(jet1, jet2) -> Float64

Function to calculate the distance in the y-ϕ plane between two jets `jet1` and `jet2`

# Arguments
- `jet1::PseudoJet`: The first jet.
- `jet2::PseudoJet`: The second jet.

# Returns
- `Float64`: The Euclidean distance in the y-ϕ plane for the two jets.
"""
function deltaR(jet1::PseudoJet, jet2::PseudoJet)
    y1, phi1 = rapidity(jet1), phi(jet1)
    y2, phi2 = rapidity(jet2), phi(jet2)

    d_y = y1 - y2
    d_phi = phi1 - phi2
    d_phi = abs(d_phi) > π ? 2π - abs(d_phi) : d_phi

    return sqrt(d_y^2 + d_phi^2)
end

"""
    recluster(jet, clusterseq; R = rad, algorithm = mtd) -> ClusterSequence

Reclusters the constituents of a given jet `jet` with a different clustering method `mtd` and different jet radius `rad`.

# Arguments
- `jet::PseudoJet`: The jet whose constituents are to be reclustered.
- `clusterseq::ClusterSequence`: The cluster sequence from which the original jet is obtained.
- `rad::Float64`: The new jet radius, by default set to 1.0
- `mtd::JetAlgorithm.Algorithm`: The new clustering method, by default set to JetAlgorithm.CA

# Returns
- `ClusterSequence`: The new cluster sequence.
"""
function recluster(jet::PseudoJet, clusterseq::ClusterSequence; R = 1.0,
                   algorithm = JetAlgorithm.CA)
    cons = constituents(jet, clusterseq)
    new_clusterseq = jet_reconstruct(cons; p = nothing, R = R, algorithm = algorithm,
                                     strategy = RecoStrategy.Best)

    return new_clusterseq
end

"""
    struct MassDropTagger

Used for tagging jets that undergo mass drop, a common technique in jet substructure.

# Fields:
- `mu::Float64`: Maximum allowed mass ratio for a jet to pass tagging
- `y::Float64`: Minimum kT distance threshold for parent separation
"""
struct MassDropTagger
    mu::Float64
    y::Float64
end

"""
    struct SoftDropTagger

Applies a soft-drop condition on jets, trimming away soft, wide-angle radiation.

# Fields:
- `zcut::Float64`: Minimum allowed energy fraction for subjets
- `b::Float64`: Angular exponent controlling soft radiation suppression
- `cluster_rad::Float64`: The new radius that will be used to recluster the components of the jet, by default set to 1.0
"""
struct SoftDropTagger
    zcut::Float64
    b::Float64
    cluster_rad::Float64
end

SoftDropTagger(z::Float64, b::Float64) = SoftDropTagger(z, b, 1.0)

"""
    struct JetFilter

Filters jets based on radius and number of hardest subjets, reducing contamination.

# Fields:
- `filter_radius::Float64`: Radius parameter to recluster subjets
- `num_hardest_jets::Int64`: Number of hardest jets to retain in the filtered result
"""
struct JetFilter
    filter_radius::Float64
    num_hardest_jets::Int64
end

"""
    struct JetTrim

Trims soft, large-angle components from jets based on fraction and radius.

# Fields:
- `trim_radius::Float64`: Radius used for reclustering in trimming
- `trim_fraction::Float64`: Minimum momentum fraction for retained subjets
- `recluster_method::JetAlgorithm.Algorithm`: Method identifier for reclustering
"""
struct JetTrim
    trim_radius::Float64
    trim_fraction::Float64
    recluster_method::JetAlgorithm.Algorithm
end

"""
    mass_drop(jet, clusterseq, tag) -> PseudoJet

Identifies subjets in a jet that pass the mass drop tagging condition.
The method stops at the first jet satisfying the mass and distance thresholds.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance representing the jet to tag
- `clusterseq::ClusterSequence`: ClusterSequence with jet clustering history
- `tag::MassDropTagger`: MassDropTagger instance providing mass drop parameters

# Returns:
- `PseudoJet`: The jet (or subjet) satisfying the mass drop conditions, if tagging is successful, otherwise a zero-momentum PseudoJet
"""
function mass_drop(jet::PseudoJet, clusterseq::ClusterSequence, tag::MassDropTagger)
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

            if m2(parent1) < m2(jet) * tag.mu^2 &&
               (min(pt1, pt2) * deltaR(parent1, parent2))^2 > tag.y * m2(jet)
                return jet
            else
                jet = parent1
            end

        else
            return PseudoJet(0.0, 0.0, 0.0, 0.0)
        end
    end
end

"""
    soft_drop(jet, clusterseq, tag) -> PseudoJet

Applies soft-drop grooming to remove soft, wide-angle radiation from jets.
This function reclusters the jet and iteratively checks the soft-drop condition on subjets.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance to groom
- `clusterseq::ClusterSequence`: ClusterSequence containing jet history
- `tag::SoftDropTagger`: SoftDropTagger instance with soft-drop parameters

# Returns:
- `PseudoJet`: Groomed jet or zero-momentum PseudoJet if grooming fails
"""
function soft_drop(jet::PseudoJet, clusterseq::ClusterSequence,
                   tag::SoftDropTagger)
    rad = tag.cluster_rad
    new_clusterseq = recluster(jet, clusterseq; R = rad, algorithm = JetAlgorithm.CA)
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
               tag.zcut * (deltaR(parent1, parent2) / rad)^tag.b
                return new_jet
            else
                new_jet = parent1
            end

        else
            return nothing
        end
    end
end

"""
    jet_filtering(jet, clusterseq, filter) -> PseudoJet

Filters a jet to retain only the hardest subjets based on a specified radius and number.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance representing the jet to filter
- `clusterseq::ClusterSequence`: ClusterSequence containing jet history
- `filter::JetFilter`: Filter instance specifying radius and number of subjets

# Returns:
- `PseudoJet`: Filtered jet composed of the hardest subjets
"""
function jet_filtering(jet::PseudoJet, clusterseq::ClusterSequence, filter::JetFilter)
    rad = filter.filter_radius
    new_clusterseq = recluster(jet, clusterseq; R = rad, algorithm = JetAlgorithm.CA)
    reclustered = sort!(inclusive_jets(new_clusterseq; T = PseudoJet), by = pt2, rev = true)

    n = length(reclustered) <= filter.num_hardest_jets ? length(reclustered) :
        filter.num_hardest_jets
    hard = reclustered[1:n]

    filtered = foldl(+, hard)

    filtered
end

"""
    jet_trimming(jet, clusterseq, trim) -> PseudoJet

Trims a jet by removing subjets with transverse momentum below a specified fraction.

# Arguments:
- `jet::PseudoJet`: PseudoJet instance representing the jet to trim
- `clusterseq::ClusterSequence`: ClusterSequence containing jet history
- `trim::JetTrim`: Trim instance specifying trimming parameters

# Returns:
- `PseudoJet`: Trimmed jet composed of retained subjets
"""
function jet_trimming(jet::PseudoJet, clusterseq::ClusterSequence, trim::JetTrim)
    rad = trim.trim_radius
    alg = trim.recluster_method
    frac2 = trim.trim_fraction^2

    new_clusterseq = recluster(jet, clusterseq; R = rad, algorithm = alg)
    reclustered = sort!(inclusive_jets(new_clusterseq; T = PseudoJet), by = pt2, rev = true)

    hard = Vector{PseudoJet}(undef, 0)
    for item in reclustered
        if pt2(item) >= frac2 * pt2(jet)
            push!(hard, item)
        end
    end
    trimmed = length(hard) != 0 ? foldl(+, hard) : PseudoJet(0.0, 0.0, 0.0, 0.0)

    trimmed
end
