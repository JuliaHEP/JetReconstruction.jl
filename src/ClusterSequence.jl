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
    struct HistoryElement{T <: Real}

A struct holding a record of jet mergers and finalisations. Floating
point fields are parameterised on T.

# Fields
- `parent1::Int`: Index in history where first parent of this jet was created
  (`NonexistentParent` if this jet is an original particle).
- `parent2::Int`: Index in history where second parent of this jet was created
  (`NonexistentParent` if this jet is an original particle; `BeamJet` if this
  history entry just labels the fact that the jet has recombined with the beam).
- `child::Int`: Index in history where the current jet is recombined with another jet
  to form its child. It is `Invalid` if this jet does not further recombine.
- `jetp_index::Int`: Index in the jets vector where we will find the `FourMomentum` object
  corresponding to this jet (i.e. the jet created at this entry of the history).
  NB: if this element of the history corresponds to a beam recombination, then
  `jetp_index=Invalid`.
- `dij::T`: The distance corresponding to the recombination at this stage of the
  clustering.
- `max_dij_so_far::T`: The largest recombination distance seen so far in the
  clustering history.
"""
struct HistoryElement{T <: Real}
    parent1::Int
    parent2::Int
    child::Int
    jetp_index::Int
    dij::T
    max_dij_so_far::T
end

import Base.eltype
"""
    eltype(::HistoryElement{T})

Return the parameterised type, `T` of the history sequence.
"""
eltype(::HistoryElement{T}) where {T} = T

"""
    HistoryElement{T}(jetp_index) where {T <: Real}

Constructs a `HistoryElement{T}` object with the given `jetp_index`, used for
initialising the history with original particles.

# Arguments
- `jetp_index`: The index of the jetp.

# Returns
A `HistoryElement{T}` object.
"""
function HistoryElement{T}(jetp_index) where {T <: Real}
    HistoryElement(NonexistentParent, NonexistentParent, Invalid,
                   jetp_index, zero(T), zero(T))
end

"""
    initial_history(particles::AbstractVector{A}) where {T <: Real, J <: FourMomentum{T},
                                              A <: AbstractVector{J}}

Create an initial history for the given particles.

# Arguments
- `particles`: The initial vector of stable particles.

# Returns
- `history`: An array of `HistoryElement{T}` objects.
- `Qtot`: The total energy in the event.
"""
function initial_history(particles::A) where {T <: Real, J <: FourMomentum{T},
                                              A <: AbstractVector{J}}
    # reserve sufficient space for everything
    history = Vector{HistoryElement{T}}(undef, length(particles))
    sizehint!(history, 2 * length(particles))

    Qtot::T = zero(T)

    for i in eachindex(particles)
        history[i] = HistoryElement{T}(i)

        # get cross-referencing right from the Jets
        # particles[i]._cluster_hist_index = i
        @assert cluster_hist_index(particles[i])==i "Cluster history index should match jet's index in the input vector. Expected $(i), got $(cluster_hist_index(particles[i]))"

        # determine the total energy in the event
        Qtot += particles[i].E
    end
    history, Qtot
end

"""
    struct ClusterSequence{T <: Real, J <: FourMomentum{T}}

A struct holding the full history of a jet clustering sequence, including the
final jets. Parameterised on the jet type `J`, including it's numerical
parameter `T`.

# Fields
- `algorithm::JetAlgorithm.Algorithm`: The algorithm used for clustering.
- `power::Float64`: The power value used for the clustering algorithm (note that
  this value is always stored as a `Float64` to be type stable).
- `R::Float64`: The R parameter used for the clustering algorithm.
- `strategy::RecoStrategy.Strategy`: The strategy used for clustering.
- `jets::Vector{J}`: The actual jets in the cluster sequence, which are of type
  `J <: FourMomentum{T}`.
- `n_initial_jets::Int`: The initial number of particles used for exclusive
  jets.
- `history::Vector{HistoryElement{T}}`: The branching history of the cluster
  sequence. Each stage in the history indicates where to look in the `jets` vector
  to get the physical jet.
- `Qtot::T`: The total energy of the event.
"""
struct ClusterSequence{T <: Real, J <: FourMomentum{T}}
    algorithm::JetAlgorithm.Algorithm
    power::T
    R::T
    strategy::RecoStrategy.Strategy
    jets::Vector{J}
    n_initial_jets::Int
    history::Vector{HistoryElement{T}}
    Qtot::T
end

import Base.eltype
"""
    eltype(::ClusterSequence{T})

Return the parameterised type, `J` of the history sequence. Note that this is not
the numeric type of the jets, that requires `eltype(eltype(::ClusterSequence))`.
"""
eltype(::ClusterSequence{T, J}) where {T, J} = J

"""
    ClusterSequence(algorithm, p, R, strategy, jets::Vector{J}, history, Qtot) where {T <: Real, J <: FourMomentum{T}}

Construct a `ClusterSequence{T, J}` object.

# Arguments
- `algorithm::JetAlgorithm.Algorithm`: The algorithm used for clustering.
- `p`: The power value used for the clustering algorithm.
- `R`: The R parameter used for the clustering algorithm.
- `strategy::RecoStrategy.Strategy`: The strategy used for clustering.
- `jets::Vector{J}`: The jets in the cluster sequence, where `J <: FourMomentum{T}`.
- `history::Vector{HistoryElement{T}}`: The branching history of the cluster
  sequence.
- `Qtot::T`: The total energy of the event.
"""
function ClusterSequence(algorithm::JetAlgorithm.Algorithm, p, R,
                         strategy::RecoStrategy.Strategy, jets::Vector{J}, history,
                         Qtot) where {T <: Real, J <: FourMomentum{T}}
    ClusterSequence{T, J}(algorithm, p, R, strategy, jets, length(jets), history, Qtot)
end

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
function add_step_to_history!(clusterseq::ClusterSequence{T, J}, parent1, parent2,
                              jetp_index,
                              dij) where {T <: Real, J <: FourMomentum{T}}
    max_dij_so_far = max(dij, clusterseq.history[end].max_dij_so_far)
    push!(clusterseq.history,
          HistoryElement{T}(parent1, parent2, Invalid,
                            jetp_index, dij, max_dij_so_far))

    local_step = length(clusterseq.history)

    # println("Adding step $local_step: parent1=$parent1, parent2=$parent2, jetp_index=$jetp_index, dij=$dij")

    # Sanity check: make sure the particles have not already been recombined
    #
    # Note that good practice would make this an assert (since this is
    # a serious internal issue). However, we decided to throw an
    # InternalError so that the end user can decide to catch it and
    # retry the clustering with a different strategy.
    @assert parent1 >= 1
    if clusterseq.history[parent1].child != Invalid
        error("Internal error. Trying to recombine an object that has previously been recombined. Parent " *
              string(parent1) * "'s child index " *
              string(clusterseq.history[parent1].child) * ". Parent jet index: " *
              string(clusterseq.history[parent1].jetp_index) * ".")
    end

    hist_elem = clusterseq.history[parent1]
    clusterseq.history[parent1] = @set hist_elem.child = local_step

    if parent2 >= 1
        clusterseq.history[parent2].child == Invalid ||
            error("Internal error. Trying to recombine an object that has previously been recombined.  Parent " *
                  string(parent2) * "'s child index " *
                  string(clusterseq.history[parent1].child) * ". Parent jet index: " *
                  string(clusterseq.history[parent2].jetp_index) * ".")
        hist_elem = clusterseq.history[parent2]
        clusterseq.history[parent2] = @set hist_elem.child = local_step
    end
end

"""
    inclusive_jets(clusterseq::ClusterSequence{T, J}, ::Type{R} = LorentzVector{Float64}; ptmin = 0.0) where {T, J, R}

Return all inclusive jets of a `ClusterSequence` with pt > ptmin.

# Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the
  clustering history and jets.
- `::Type{R} = LorentzVector{Float64}`: The return type used for the selected jets.
- `ptmin::Real = 0.0`: The minimum transverse momentum (pt) threshold for the
  inclusive jets.

# Returns
An array of `R` objects representing the inclusive jets.

# Description
This function computes the inclusive jets from a given `ClusterSequence` object.
It iterates over the clustering history and checks the transverse momentum of
each jet that merged with the beam. If the transverse momentum is greater than 
or equal to `ptmin`, the jet is added to the array of inclusive jets.

Valid return types are `LorentzVector`, `LorentzVectorCyl`, or the jet type of the
input `clusterseq` (`J` - usually either `PseudoJet` or `EEJet`).

# Example
```julia
inclusive_jets(clusterseq; ptmin = 10.0)
```
"""
function inclusive_jets(clusterseq::ClusterSequence{T, J},
                        ::Type{R} = LorentzVector{T};
                        ptmin = 0.0) where {T <: Real, J <: FourMomentum{T}, R}
    pt2min = ptmin * ptmin
    RetType = concretize_return_type(R, T)
    jets_local = RetType[]
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
            if J == RetType
                push!(jets_local, jet)
            elseif RetType <: LorentzVectorCyl
                push!(jets_local, lorentzvector_cyl(jet))
            elseif RetType <: LorentzVector
                push!(jets_local, lorentzvector(jet))
            else
                error("Unsupported return type $RetType for inclusive jets (cf. ClusterSequence jet type $J)")
            end
        end
    end
    jets_local
end

"""
    exclusive_jets(clusterseq::ClusterSequence{T, J}, ::Type{R} = LorentzVector{Float64}; dcut = nothing, njets = nothing) where {T, J, R}

Return all exclusive jets of a `ClusterSequence`, with either a specific number of
jets or a cut on the maximum distance parameter.

# Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the
  clustering history and jets.
- `::Type{R} = LorentzVector{Float64}`: The return type used for the selected
  jets.
- `dcut::Union{Nothing, Real}`: The distance parameter used to define the
  exclusive jets. If `dcut` is provided, the number of exclusive jets will be
  calculated based on this parameter.
- `njets::Union{Nothing, Integer}`: The number of exclusive jets to be
  calculated. If `njets` is provided, the distance parameter `dcut` will be
  calculated based on this number.

**Note**: Either `dcut` or `njets` must be provided (but not both).

# Returns
- An array of `R` objects representing the exclusive jets.

Valid return types are `LorentzVector`, `LorentzVectorCyl`, or the jet type of the
input `clusterseq` (`J` - usually either `PseudoJet` or `EEJet`).

# Exceptions
- `ArgumentError`: If neither `dcut` nor `njets` is provided.
- `ArgumentError`: If the algorithm used in the `ClusterSequence` object is not
  suitable for exclusive jets.
- `ErrorException`: If the cluster sequence is incomplete and exclusive jets are
  unavailable.

# Examples
```julia
exclusive_jets(clusterseq, dcut = 20.0)
exclusive_jets(clusterseq, PseudoJet, njets = 3)
```
"""
function exclusive_jets(clusterseq::ClusterSequence{T, J},
                        ::Type{R} = LorentzVector{Float64};
                        dcut = nothing,
                        njets = nothing,) where {T <: Real, J <: FourMomentum{T}, R}
    if isnothing(dcut) && isnothing(njets)
        throw(ArgumentError("Must pass either a dcut or an njets value"))
    end

    if !isnothing(dcut)
        njets = n_exclusive_jets(clusterseq, dcut = dcut)
    end

    # Check that an algorithm was used that makes sense for exclusive jets
    if (clusterseq.algorithm ∈ (JetAlgorithm.GenKt, JetAlgorithm.EEKt)) &&
       clusterseq.power < 0
        throw(ArgumentError("Algorithm $(clusterseq.algorithm) requires power >= 0 for exclusive jets (power=$(clusterseq.power))"))
    elseif clusterseq.algorithm ∉
           (JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.Durham, JetAlgorithm.GenKt,
            JetAlgorithm.EEKt, JetAlgorithm.Valencia)
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

    RetType = concretize_return_type(R, T)
    excl_jets = RetType[]
    for j in stop_point:length(clusterseq.history)
        @debug "Search $j ($(clusterseq.history[j].parent1) + $(clusterseq.history[j].parent2))"
        for parent in (clusterseq.history[j].parent1, clusterseq.history[j].parent2)
            if (parent < stop_point && parent > 0)
                @debug "Added exclusive jet index $(clusterseq.history[parent].jetp_index)"
                jet = clusterseq.jets[clusterseq.history[parent].jetp_index]
                if J == RetType
                    push!(excl_jets, jet)
                elseif RetType <: LorentzVectorCyl
                    push!(excl_jets, lorentzvector_cyl(jet))
                elseif RetType <: LorentzVector
                    push!(excl_jets, lorentzvector(jet))
                else
                    error("Unsupported return type $RetType for exclusive jets (cf. ClusterSequence jet type $J)")
                end
            end
        end
    end
    excl_jets
end

"""
    n_exclusive_jets(clusterseq::ClusterSequence, dcut::Real)

Return the number of exclusive jets of a `ClusterSequence` that are above a certain `dcut` value.

# Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the clustering history.
- `dcut::Real`: The maximum value for the distance parameter in the reconstruction.

# Returns
The number of exclusive jets in the `ClusterSequence` object.

# Example
```julia
n_exclusive_jets(clusterseq, dcut = 20.0)
```
"""
function n_exclusive_jets(clusterseq::ClusterSequence; dcut::AbstractFloat)
    # Check that an algorithm was used that makes sense for exclusive jets
    if (clusterseq.algorithm ∈ (JetAlgorithm.GenKt, JetAlgorithm.EEKt)) &&
       clusterseq.power < 0
        throw(ArgumentError("Algorithm $(clusterseq.algorithm) requires power >= 0 for exclusive jets(power=$(clusterseq.power))"))
    elseif clusterseq.algorithm ∉
           (JetAlgorithm.CA, JetAlgorithm.Kt, JetAlgorithm.Durham, JetAlgorithm.GenKt,
            JetAlgorithm.EEKt, JetAlgorithm.Valencia)
        throw(ArgumentError("Algorithm used is not suitable for exclusive jets ($(clusterseq.algorithm))"))
    end

    # Locate the point where clustering would have stopped (i.e. the
    # first time max_dij_so_far > dcut)
    i_dcut = length(clusterseq.history)
    for i_history in length(clusterseq.history):-1:1
        @debug "Examining $i_history, max_dij=$(clusterseq.history[i_history].max_dij_so_far)"
        if clusterseq.history[i_history].max_dij_so_far <= dcut
            i_dcut = i_history
            break
        end
    end

    # The number of jets is then given by this formula
    length(clusterseq.history) - i_dcut
end

"""
    get_all_ancestors(idx, cs::ClusterSequence)

Recursively finds all ancestors of a given index in a `ClusterSequence` object.

# Arguments
- `idx`: The index of the jet for which to find ancestors.
- `cs`: The `ClusterSequence` object containing the jet history.

# Returns
An array of indices representing the ancestors of the given jet.
"""
function get_all_ancestors(idx, cs::ClusterSequence)
    if cs.history[idx].parent1 == JetReconstruction.NonexistentParent
        return [cs.history[idx].jetp_index]
    else
        branch1 = get_all_ancestors(cs.history[idx].parent1, cs)
        cs.history[idx].parent2 == JetReconstruction.BeamJet && return branch1
        branch2 = get_all_ancestors(cs.history[idx].parent2, cs)
        return [branch1; branch2]
    end
end

"""
    merge_steps(clusterseq::ClusterSequence)

Compute the number of jet-jet merge steps in a cluster sequence. This is useful
to give the number of meaningful recombination steps in a jet reconstruction
sequence (beam merge steps are not counted).

# Arguments
- `clusterseq::ClusterSequence`: The cluster sequence object.

# Returns
- `merge_steps::Int`: The number of merge steps.
"""
function merge_steps(clusterseq::ClusterSequence)
    merge_steps = 0
    for step in clusterseq.history[(clusterseq.n_initial_jets + 1):end]
        step.parent2 == BeamJet && continue
        merge_steps += 1
    end
    merge_steps
end

"""
    jet_ranks(clusterseq::ClusterSequence; compare_fn = JetReconstruction.pt)

Compute the ranks of jets in a given `ClusterSequence` object based on a
specified comparison function.

## Arguments
- `clusterseq::ClusterSequence`: The `ClusterSequence` object containing the
  jets to rank.
- `compare_fn = JetReconstruction.pt`: The comparison function used to determine
  the order of the jets. Defaults to `JetReconstruction.pt`, which compares jets
  based on their transverse momentum.

## Returns
A dictionary mapping each jet index to its rank.

## Note
This is a utility function that can be used to rank initial clusters based on a specified
jet property. It can be used to assign a consistent "rank" to each reconstructed jet in
the cluster sequence, which is useful for stable plotting of jet outputs.
"""
function jet_ranks(clusterseq::ClusterSequence; compare_fn = JetReconstruction.pt)
    initial_jet_list = collect(1:(clusterseq.n_initial_jets))
    sort!(initial_jet_list, by = i -> compare_fn(clusterseq.jets[i]), rev = true)
    jet_ranks = Dict{Int, Int}()
    for (rank, jetp_index) in enumerate(initial_jet_list)
        jet_ranks[jetp_index] = rank
    end
    jet_ranks
end

"""
    struct JetWithAncestors{J <: FourMomentum}

A struct representing a jet with its origin ancestors.

# Fields
- `self::J`: The jet object itself.
- `jetp_index::Int`: The index of the jet in the corresponding cluster sequence.
- `ancestors::Set{Int}`: A set of indices representing the `jetp_index`es of
  ancestors of the jet (in the cluster sequence).
- `jet_rank::Int`: The rank of the jet based on a comparison of all of the jet's
  ancestors.

# Note
This structure needs its associated cluster sequence origin to be useful.
"""
struct JetWithAncestors{J <: FourMomentum}
    self::J
    jetp_index::Int
    ancestors::Set{Int}
    jet_rank::Int
end

"""
    reco_state(cs::ClusterSequence{T, J}, ranks; iteration = 0, ignore_beam_merge = true) where {T, J}

This function returns the reconstruction state of a `ClusterSequence` object
based on a given iteration number in the reconstruction.

# Arguments
- `cs::ClusterSequence`: The `ClusterSequence` object to update.
- `ranks`: The ranks of the original clusters, that are inherited by pseudojets
 during the reconstruction process.
- `iteration = 0`: The iteration number to consider for updating the
  reconstruction state (0 represents the initial state).
- `ignore_beam_merge = true`: Ignore beam merging steps in the reconstruction
  (which produce no change in status).

# Returns
A dictionary representing a snapshot of the reconstruction state, mapping 
jet indices to `JetWithAncestors` objects.

# Details
The function starts by initializing the reconstruction state with the initial
particles. Then, it walks over the iteration sequence and updates the
reconstruction state based on the history of recombination and finalization/beam
merger steps.
"""
function reco_state(cs::ClusterSequence{T, J}, ranks; iteration = 0,
                    ignore_beam_merge = true) where {T <: Real, J <: FourMomentum{T}}
    # Get the initial particles
    reco_state = Dict{Int, JetWithAncestors{J}}()
    for jet_index in 1:(cs.n_initial_jets)
        reco_state[jet_index] = JetWithAncestors{J}(cs.jets[cs.history[jet_index].jetp_index],
                                                    jet_index, Set{Int}([]),
                                                    ranks[jet_index])
    end
    # Now update the reconstruction state by walking over the iteration sequence
    iterations_done = 0
    for h_step in (cs.n_initial_jets + 1):length(cs.history)
        h_entry = cs.history[h_step]
        if h_entry.parent2 > 0
            # This is a recombination
            iterations_done += 1
            # We store all of the original particle ancestors (but only the
            # original ones, not intermediate merges)
            my_ancestors = union(reco_state[cs.history[h_entry.parent1].jetp_index].ancestors,
                                 reco_state[cs.history[h_entry.parent2].jetp_index].ancestors)
            cs.history[h_entry.parent1].parent1 == JetReconstruction.NonexistentParent &&
                push!(my_ancestors, h_entry.parent1)
            cs.history[h_entry.parent2].parent1 == JetReconstruction.NonexistentParent &&
                push!(my_ancestors, h_entry.parent2)

            # Now find the ancestor with the highest p_T value
            pt_rank = cs.n_initial_jets
            for ancestor in my_ancestors
                (ranks[ancestor] < pt_rank) && (pt_rank = ranks[ancestor])
            end

            reco_state[h_entry.jetp_index] = JetWithAncestors{T}(cs.jets[h_entry.jetp_index],
                                                                 h_entry.jetp_index,
                                                                 my_ancestors, pt_rank)
            delete!(reco_state, cs.history[h_entry.parent1].jetp_index)
            delete!(reco_state, cs.history[h_entry.parent2].jetp_index)
        else
            # This is a finalisation / beam merger, so here we do nothing
            if ignore_beam_merge != true
                iterations_done += 1
            end
        end
        # Abort when we have done the required number of iterations
        if iterations_done == iteration
            break
        end
    end
    reco_state
end

"""
    constituents(jet::J, cs::ClusterSequence{T, J}) where {T, J}

Get a copy of the constituents of a given jet in a cluster sequence.

# Arguments
- `jet::J`: The jet for which to retrieve the constituents.
- `cs::ClusterSequence{T, J}`: The cluster sequence object.

# Returns
An array of jet objects (of type `J`) copied from the constituents of the given 
jet, with reset cluster history indexes.
"""
function constituents(jet::J,
                      cs::ClusterSequence{T, J}) where {T <: Real, J <: FourMomentum{T}}
    constituent_idxs = constituent_indexes(jet, cs)
    constituents = Vector{J}()
    sizehint!(constituents, length(constituent_idxs))
    new_index = 1
    for idx in constituent_idxs
        push!(constituents,
              J(px(cs.jets[idx]), py(cs.jets[idx]), pz(cs.jets[idx]),
                energy(cs.jets[idx]); cluster_hist_index = new_index))
        new_index += 1
    end
    constituents
end

"""
    constituent_indexes(jet::J, cs::ClusterSequence{T, J}) where {T, J}

Return the indexes of the original particles which are the constituents of the
given jet.

# Arguments
- `jet::J`: The jet for which to retrieve the constituents.
- `cs::ClusterSequence{T, J}`: The cluster sequence object.

# Returns
A vector of indices representing the original constituents of the given jet.
"""
function constituent_indexes(jet::J,
                             cs::ClusterSequence{T, J}) where {T <: Real,
                                                               J <: FourMomentum{T}}
    get_all_ancestors(jet._cluster_hist_index, cs)
end

"""
    parent_jets(jet::J, cs::ClusterSequence{T, J}) where {T, J}

Find the parent jets of a given jet in a cluster sequence.

# Arguments
- `jet::J`: The jet for which to find the parent jets.
- `cs::ClusterSequence{T, J}`: The cluster sequence object.

# Returns
A tuple of two elements, each of which is either the parent jet object (of type `J`) 
or `nothing` (if the jet has no parent).
"""
function parent_jets(jet::J,
                     cs::ClusterSequence{T, J})::Tuple{Union{Nothing, J},
                                                       Union{Nothing, J}} where {T <:
                                                                                 Real,
                                                                                 J <:
                                                                                 FourMomentum{T}
                                                                                 }
    hist_idx = jet._cluster_hist_index
    jet_history = cs.history[hist_idx]

    parent1_idx, parent2_idx = jet_history.parent1, jet_history.parent2

    parent1_jet = parent1_idx > 0 ? cs.jets[cs.history[parent1_idx].jetp_index] : nothing
    parent2_jet = parent2_idx > 0 ? cs.jets[cs.history[parent2_idx].jetp_index] : nothing

    return parent1_jet, parent2_jet
end
