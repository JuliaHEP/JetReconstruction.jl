"""
    has_parents(p, history) -> (boolean, Int64, Int64)

Checks if the jet `p` is a child of two other jets, after clustering 

# Arguments
- `p`: The jet to check.
- `history`: The vector of history element obtained from the cluster sequence after clustering.

# Returns
- (boolean, Int64, Int64): true or false depending on if the jet has a parent or not. If the jet has a parent, returns the indices of the parent jets in the history element. Otherwise, returns -2 (NonexistentParent).
"""
function has_parents(p::PseudoJet, history::Vector{HistoryElement})
    
    N = p._cluster_hist_index
    p1 = history[N].parent1
    p2 = history[N].parent2

    p1 == p2 == NonexistentParent ? result = false : result = true
    return (result, p1, p2)

end

"""
    delta_R(jet1, jet2) -> Float64

Function to calculate the distance in the y-ϕ plane between two jets `jet1` and `jet2`

# Arguments
- `jet1`: The first jet.
- `jet2`: The second jet.

# Returns
- Float64: The Euclidean distance in the y-ϕ plane for the two jets.
"""
function delta_R(jet1::PseudoJet, jet2::PseudoJet)
    
    eta1, phi1 = rapidity(jet1), phi(jet1)
    eta2, phi2 = rapidity(jet2), phi(jet2)

    d_eta = eta1 - eta2
    d_phi = phi1 - phi2
    d_phi = abs(d_phi) > π ? 2π - abs(d_phi) : d_phi
    
    return sqrt(d_eta^2 + d_phi^2)

end

# Function to calculate kt distance between two pseudojets
function kt_distance(jet1::PseudoJet, jet2::PseudoJet, R = 1)
    p1 = pt2(jet1)
    p2 = pt2(jet2)
    
    d_R = delta_R(jet1, jet2)
    return min(p1, p2) * (d_R^2 / R^2)
end;

"""
    recluster(jet, clusterseq, rad, mtd) -> ClusterSequence

Reclusters the constituents of a given jet `jet` with a different clustering method `mtd` and different jet radius `rad`.

# Arguments
- `jet`: The jet whose constituents are to be reclustered.
- `clusterseq`: The cluster sequence from which the original jet is obtained.
- `rad`: The new jet radius.
- `mtd`: The new clustering method.

# Returns
- ClusterSequence: The new cluster sequence.
"""
function recluster(jet::PseudoJet, clusterseq::ClusterSequence, rad = 1.0, mtd = 0)
    cons = constituents(jet, clusterseq)
    new_clusterseq = jet_reconstruct(cons; p = mtd, R = rad, strategy = RecoStrategy.N2Plain)
   
    return new_clusterseq
end

function sort_jets!(event_jet_array)
    jet_pt(jet) = pt2(jet)
    sort!(event_jet_array, by = jet_pt, rev = true)
    # println(event_jet_array)
end

function join(jets::Vector{PseudoJet})
    final = PseudoJet(0.0, 0.0, 0.0, 0.0)
    
    for jet in jets
        final = final + jet
    end

    final
end

# Defining suitable structures to be used for jet grooming and tagging
"""
    struct MassDropTagger

Used for tagging jets that undergo mass drop, a common technique in jet substructure.

Fields:
- mu: Maximum allowed mass ratio for a jet to pass tagging
- y: Minimum kT distance threshold for parent separation
"""

mutable struct MassDropTagger
    mu::Float64
    y::Float64
end

"""
    struct SoftDropTagger

Applies a soft-drop condition on jets, trimming away soft, wide-angle radiation.

Fields:
- zcut: Minimum allowed energy fraction for subjets
- b: Angular exponent controlling soft radiation suppression
"""

mutable struct SoftDropTagger
    zcut::Float64
    b::Float64
end

"""
    struct Filter

Filters jets based on radius and number of hardest subjets, reducing contamination.

Fields:
- filterRadius: Radius parameter to recluster subjets
- numHardestJets: Number of hardest jets to retain in the filtered result
"""

mutable struct Filter
    filterRadius::Float64 
    numHardestJets::Int
end

"""
    struct Trim

Trims soft, large-angle components from jets based on fraction and radius.

Fields:
- trimRadius: Radius used for reclustering in trimming
- trimFraction: Minimum momentum fraction for retained subjets
- reclusterMethod: Method identifier for reclustering
"""

mutable struct Trim
    trimRadius::Float64
    trimFraction::Float64
    reclusterMethod::Int
end

"""
    massDrop(jet, clusterseq, tag) -> PseudoJet

Identifies subjets in a jet that pass the mass drop tagging condition.
The method stops at the first jet satisfying the mass and distance thresholds.

Arguments:
- jet: PseudoJet instance representing the jet to tag
- clusterseq: ClusterSequence with jet clustering history
- tag: MassDropTagger instance providing mass drop parameters

Returns:
- PseudoJet: The jet (or subjet) satisfying the mass drop conditions, if tagging is successful, otherwise a zero-momentum PseudoJet
"""

function massDrop(jet::PseudoJet, clusterseq::ClusterSequence, tag::MassDropTagger)
    allJets = clusterseq.jets
    hist = clusterseq.history

    while(true)
        had_parents, p1, p2 = has_parents(jet, hist)

        if (had_parents)
            parent1 = allJets[hist[p1].jetp_index]
            parent2 = allJets[hist[p2].jetp_index]

            if (m2(parent1) < m2(parent2))
                p1, p2 = p2, p1
                parent1, parent2 = parent2, parent1
            end

            if ((m2(parent1) < m2(jet)*tag.mu^2) && (kt_distance(parent1, parent2) > tag.y*m2(jet)))
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
    softDrop(jet, clusterseq, rad, tag) -> PseudoJet

Applies soft-drop grooming to remove soft, wide-angle radiation from jets.
This function reclusters the jet and iteratively checks the soft-drop condition on subjets.

Arguments:
- jet: PseudoJet instance to groom
- clusterseq: ClusterSequence containing jet history
- rad: Radius parameter for reclustering
- tag: SoftDropTagger instance with soft-drop parameters

Returns:
- PseudoJet: Groomed jet or zero-momentum PseudoJet if grooming fails
"""

function softDrop(jet::PseudoJet, clusterseq::ClusterSequence, rad::Float64, tag::SoftDropTagger)
    new_clusterseq = recluster(jet, clusterseq, rad, 0)
    new_jet = inclusive_jets(new_clusterseq; T = PseudoJet)[1]

    allJets = new_clusterseq.jets
    hist = new_clusterseq.history
    
    while(true)
        had_parents, p1, p2 = has_parents(new_jet, hist)

        if (had_parents)
            parent1 = allJets[hist[p1].jetp_index]
            parent2 = allJets[hist[p2].jetp_index]

            pti = pt2(parent1)^0.5
            ptj = pt2(parent2)^0.5

            if (m2(parent1) < m2(parent2))
                p1, p2 = p2, p1
                parent1, parent2 = parent2, parent1
            end

            if (min(pti, ptj)/(pti + ptj) > tag.zcut * (delta_R(parent1, parent2)/rad) ^ tag.b)
                return new_jet
            else
                new_jet = parent1
            end
    
        else 
            return PseudoJet(0.0, 0.0, 0.0, 0.0)
        end

    end

end

"""
    jetFilter(jet, clusterseq, filter) -> PseudoJet

Filters a jet to retain only the hardest subjets based on a specified radius and number.

Arguments:
- jet: PseudoJet instance representing the jet to filter
- clusterseq: ClusterSequence containing jet history
- filter: Filter instance specifying radius and number of subjets

Returns:
- PseudoJet: Filtered jet composed of the hardest subjets

"""
function jetFilter(jet::PseudoJet, clusterseq::ClusterSequence, filter::Filter)
    rad = filter.filterRadius;
    new_clusterseq = recluster(jet, clusterseq, rad, 0)
    reclustered = sort_jets!(inclusive_jets(new_clusterseq; T = PseudoJet))

    n = length(reclustered) <= filter.numHardestJets ? length(reclustered) : filter.numHardestJets 
    hard = reclustered[1:n]
    
    filtered = join(hard)
    
    filtered    
end

"""
    jetTrim(jet, clusterseq, trim) -> PseudoJet

Trims a jet by removing subjets with transverse momentum below a specified fraction.

Arguments:
- jet: PseudoJet instance representing the jet to trim
- clusterseq: ClusterSequence containing jet history
- trim: Trim instance specifying trimming parameters

Returns:
- PseudoJet: Trimmed jet composed of retained subjets
"""

function jetTrim(jet::PseudoJet, clusterseq::ClusterSequence, trim::Trim)
    rad = trim.trimRadius;
    mtd = trim.reclusterMethod
    frac2 = trim.trimFraction ^ 2

    new_clusterseq = recluster(jet, clusterseq, rad, mtd)
    reclustered = sort_jets!(inclusive_jets(new_clusterseq; T = PseudoJet))
    
    hard = Vector{PseudoJet}(undef, 0)
    for item in reclustered
        if pt2(item) >= frac2 * pt2(jet)
            push!(hard, item)
        end
    end
    trimmed = join(hard)
    
    trimmed    
end;
