# Structure definitions for the Tiled algorithm, with linked list
# and TiledJets (a la FastJet)
#
# Original Julia implementation by Philippe Gras,
# ported to this package by Graeme Stewart

using Accessors
using Logging

# Constant values for "magic numbers" that represent special states
# in the history. These are negative integers, but as this is 
# not runtime critical, so could consider a Union{Int,Emum}
"Invalid child for this jet, meaning it did not recombine further"
const Invalid = -3

"Original cluster so it has no parent"
const NonexistentParent = -2

"Cluster recombined with beam"
const BeamJet = -1

"""
    struct HistoryElement

A struct holding a record of jet mergers and finalisations

Fields:
- parent1: Index in history where first parent of this jet was created (NonexistentParent if this jet is an original particle)
- parent2: Index in history where second parent of this jet was created (NonexistentParent if this jet is an original particle); BeamJet if this history entry just labels the fact that the jet has recombined with the beam)
- child: Index in history where the current jet is recombined with another jet to form its child. It is Invalid if this jet does not further recombine.
- jetp_index: Index in the jets vector where we will find the PseudoJet object corresponding to this jet (i.e. the jet created at this entry of the history). NB: if this element of the history corresponds to a beam recombination, then jetp_index=Invalid.
- dij: The distance corresponding to the recombination at this stage of the clustering.
- max_dij_so_far: The largest recombination distance seen so far in the clustering history.
"""
struct HistoryElement
    parent1::Int
    parent2::Int
    child::Int
    jetp_index::Int
    dij::Float64
    max_dij_so_far::Float64
end


"""
    HistoryElement(jetp_index)

Constructs a `HistoryElement` object with the given `jetp_index`, used for initialising the history with original particles.

# Arguments
- `jetp_index`: The index of the jetp.

# Returns
A `HistoryElement` object.

"""
HistoryElement(jetp_index) = HistoryElement(NonexistentParent, NonexistentParent, Invalid, jetp_index, 0.0, 0.0)


"""
    struct ClusterSequence

A struct holding the full history of a jet clustering sequence, including the final jets.

# Fields
- `algorithm::JetAlgorithm.Algorithm`: The algorithm used for clustering.
- `strategy::RecoStrategy.Strategy`: The strategy used for clustering.
- `jets::Vector{PseudoJet}`: The physical PseudoJets in the cluster sequence. Each PseudoJet corresponds to a position in the history.
- `n_initial_jets::Int`: The initial number of particles used for exclusive jets.
- `history::Vector{HistoryElement}`: The branching history of the cluster sequence. Each stage in the history indicates where to look in the jets vector to get the physical PseudoJet.
- `Qtot::Any`: The total energy of the event.
"""
struct ClusterSequence
    algorithm::JetAlgorithm.Algorithm
    strategy::RecoStrategy.Strategy
    jets::Vector{PseudoJet}
    n_initial_jets::Int
    history::Vector{HistoryElement}
    Qtot::Any
end

"""
    ClusterSequence(p::Int, strategy::RecoStrategy.Strategy, jets, history, Qtot)

Constructs a `ClusterSequence` object with a specific power value mapped to a pp reconstruction algorithm.

# Arguments
- `p::Int`: The power value for the algorithm.
- `strategy::RecoStrategy.Strategy`: The reconstruction strategy.
- `jets`: The jets to be clustered.
- `history`: The length of the jets.
- `Qtot`: The total energy of the jets.

# Returns
A `ClusterSequence` object.
"""
ClusterSequence(p::Int, strategy::RecoStrategy.Strategy, jets, history, Qtot) = begin
    if !haskey(power2algorithm, p)
        raise(ArgumentError("Unrecognised algorithm for power value p=$p"))
    end
    ClusterSequence(power2algorithm[p], strategy, jets, length(jets), history, Qtot)
end

"""
    ClusterSequence(alg::JetAlgorithm.Algorithm, strategy::RecoStrategy.Strategy, jets, history, Qtot)

Constructs a `ClusterSequence` object with a specific algorithm spcified.

# Arguments
- `alg::JetAlgorithm.Algorithm`: The algorithm.
- `strategy::RecoStrategy.Strategy`: The reconstruction strategy.
- `jets`: The jets to be clustered.
- `history`: The length of the jets.
- `Qtot`: The total energy of the jets.

# Returns
A `ClusterSequence` object.
"""
ClusterSequence(alg::JetAlgorithm.Algorithm, strategy::RecoStrategy.Strategy, jets, history, Qtot) =
    ClusterSequence(alg, strategy, jets, length(jets), history, Qtot)

"""
    add_step_to_history!(clusterseq::ClusterSequence, parent1, parent2, jetp_index, dij)

Add a new jet's history into the recombination sequence.

Arguments:
- `clusterseq::ClusterSequence`: The cluster sequence object.
- `parent1`: The index of the first parent.
- `parent2`: The index of the second parent.
- `jetp_index`: The index of the jet.
- `dij`: The dij value.

This function adds a new `HistoryElement` to the `history` vector of the
`clusterseq` object. The `HistoryElement` contains information about the
parents, child, jet index, dij value, and the maximum dij value so far. It also
updates the child index of the parent elements.

If the `parent1` or `parent2` have already been recombined, an `InternalError`
is thrown. The `jetp_index` is used to update the `_cluster_hist_index` of the
corresponding `PseudoJet` object.
"""
add_step_to_history!(clusterseq::ClusterSequence, parent1, parent2, jetp_index, dij) = begin
    max_dij_so_far = max(dij, clusterseq.history[end].max_dij_so_far)
    push!(clusterseq.history, HistoryElement(parent1, parent2, Invalid,
        jetp_index, dij, max_dij_so_far))

    local_step = length(clusterseq.history)

    # Sanity check: make sure the particles have not already been recombined
    #
    # Note that good practice would make this an assert (since this is
    # a serious internal issue). However, we decided to throw an
    # InternalError so that the end user can decide to catch it and
    # retry the clustering with a different strategy.
    @assert parent1 >= 1
    if clusterseq.history[parent1].child != Invalid
        error("Internal error. Trying to recombine an object that has previsously been recombined.")
    end

    hist_elem = clusterseq.history[parent1]
    clusterseq.history[parent1] = @set hist_elem.child = local_step

    if parent2 >= 1
        clusterseq.history[parent2].child == Invalid || error(
            "Internal error. Trying to recombine an object that has previsously been recombined.  Parent " * string(parent2) * "'s child index " * string(clusterseq.history[parent1].child) * ". Parent jet index: " *
            string(clusterseq.history[parent2].jetp_index) * ".",
        )
        hist_elem = clusterseq.history[parent2]
        clusterseq.history[parent2] = @set hist_elem.child = local_step
    end

    # Get cross-referencing right from PseudoJets
    if jetp_index != Invalid
        @assert jetp_index >= 1
        clusterseq.jets[jetp_index]._cluster_hist_index = local_step
    end
end

"""
    inclusive_jets(clusterseq::ClusterSequence, ptmin = 0.0)

Return all inclusive jets of a ClusterSequence with pt > ptmin.

# Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the clustering history and jets.
- `ptmin::Float64 = 0.0`: The minimum transverse momentum (pt) threshold for the inclusive jets.

# Returns
An array of `LorentzVectorCyl` objects representing the inclusive jets.

# Description
This function computes the inclusive jets from a given `ClusterSequence` object. It iterates over the clustering history and checks the transverse momentum of each parent jet. If the transverse momentum is greater than or equal to `ptmin`, the jet is added to the array of inclusive jets.

# Example
```julia
inclusive_jets(clusterseq, ptmin = 10.0)
```
"""
function inclusive_jets(clusterseq::ClusterSequence, ptmin = 0.0)
    pt2min = ptmin * ptmin
    jets_local = LorentzVectorCyl[]
    # sizehint!(jets_local, length(clusterseq.jets))
    # For inclusive jets with a plugin algorithm, we make no
    # assumptions about anything (relation of dij to momenta,
    # ordering of the dij, etc.)
    # for elt in Iterators.reverse(clusterseq.history)
    for elt in clusterseq.history
        elt.parent2 == BeamJet || continue
        iparent_jet = clusterseq.history[elt.parent1].jetp_index
        jet = clusterseq.jets[iparent_jet]
        if pt2(jet) >= pt2min
            @debug "Added inclusive jet index $iparent_jet"
            push!(jets_local, LorentzVectorCyl(pt(jet), rapidity(jet), phi(jet), mass(jet)))
        end
    end
    jets_local
end

"""
    exclusive_jets(clusterseq::ClusterSequence; dcut = nothing, njets = nothing)

Return all exclusive jets of a ClusterSequence, with either a specific number of jets
or a cut on the maximum distance parameter.

# Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the clustering history and jets.
- `dcut::Union{Nothing, Real}`: The distance parameter used to define the exclusive jets. If `dcut` is provided, the number of exclusive jets will be calculated based on this parameter.
- `njets::Union{Nothing, Integer}`: The number of exclusive jets to be calculated. If `njets` is provided, the distance parameter `dcut` will be calculated based on this number.

**Note**: Either `dcut` or `njets` must be provided (but not both).

# Returns
- `excl_jets::Array{LorentzVectorCyl}`: An array of `LorentzVectorCyl` objects representing the exclusive jets.

# Exceptions
- `ArgumentError`: If neither `dcut` nor `njets` is provided.
- `ArgumentError`: If the algorithm used in the `ClusterSequence` object is not suitable for exclusive jets.
- `ErrorException`: If the cluster sequence is incomplete and exclusive jets are unavailable.

# Examples
```julia
exclusive_jets(clusterseq, dcut = 20.0)
exclusive_jets(clusterseq, njets = 3)
```
"""
function exclusive_jets(clusterseq::ClusterSequence; dcut = nothing, njets = nothing)
    if isnothing(dcut) && isnothing(njets)
        throw(ArgumentError("Must pass either a dcut or an njets value"))
    end

    if !isnothing(dcut)
        njets = n_exclusive_jets(clusterseq, dcut = dcut)
    end

    # Check that an algorithm was used that makes sense for exclusive jets
    if !(clusterseq.algorithm ∈ (JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.EEKt, JetAlgorithm.Durham))
        throw(ArgumentError("Algorithm used is not suitable for exclusive jets ($(clusterseq.algorithm))"))
    end

    # njets search
    stop_point = 2 * clusterseq.n_initial_jets - njets + 1

    # Sanity check - never return more jets than initial particles
    if stop_point < clusterseq.n_initial_jets
        stop_point = clusterseq.n_initial_jets
    end

    # Sanity check - ensure that reconstruction proceeded to the end
    if 2 * clusterseq.n_initial_jets != length(clusterseq.history)
        throw(ErrorException("Cluster sequence is incomplete, exclusive jets unavailable"))
    end

    excl_jets = LorentzVectorCyl[]
    for j in stop_point:length(clusterseq.history)
        @debug "Search $j ($(clusterseq.history[j].parent1) + $(clusterseq.history[j].parent2))"
        for parent in (clusterseq.history[j].parent1, clusterseq.history[j].parent2)
            if (parent < stop_point && parent > 0)
                @debug "Added exclusive jet index $(clusterseq.history[parent].jetp_index)"
                jet = clusterseq.jets[clusterseq.history[parent].jetp_index]
                push!(excl_jets, LorentzVectorCyl(pt(jet), rapidity(jet), phi(jet), mass(jet)))
            end
        end
    end

    excl_jets
end


"""
    n_exclusive_jets(clusterseq::ClusterSequence; dcut::AbstractFloat)

Return the number of exclusive jets of a ClusterSequence that are above a certain dcut value.

# Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the clustering history.
- `dcut::AbstractFloat`: The maximum calue for the distance parameter in the reconstruction.

# Returns
The number of exclusive jets in the `ClusterSequence` object.

# Example
```julia
n_exclusive_jets(clusterseq, dcut = 20.0)
```
"""
function n_exclusive_jets(clusterseq::ClusterSequence; dcut::AbstractFloat)
    # Check that an algorithm was used that makes sense for exclusive jets
    if !(clusterseq.algorithm ∈ (JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.EEKt, JetAlgorithm.Durham))
        throw(ArgumentError("Algorithm used is not suitable for exclusive jets ($(clusterseq.algorithm))"))
    end

    # Locate the point where clustering would have stopped (i.e. the
    # first time max_dij_so_far > dcut)
    i_dcut = length(clusterseq.history)
    for i_history ∈ length(clusterseq.history):-1:1
        @debug "Examining $i_history, max_dij=$(clusterseq.history[i_history].max_dij_so_far)"
        if clusterseq.history[i_history].max_dij_so_far <= dcut
            i_dcut = i_history
            break
        end
    end

    # The number of jets is then given by this formula
    length(clusterseq.history) - i_dcut
end
