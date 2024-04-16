# Structure definitions for the Tiled algorithm, with linked list
# and TiledJets (a la FastJet)
# Original Julia implementation by Philippe Gras,
# ported to this package by Graeme Stewart

using Accessors

# One cannot use an ENUM here as it's a different type
# I don't know any good way to keep this lean and to have these necessary
# flags other than using these "magic numbers"
const Invalid=-3
const NonexistentParent=-2
const BeamJet=-1

"""A struct holding a record of jet mergers and finalisations"""
struct HistoryElement
    """Index in _history where first parent of this jet
    was created (NonexistentParent if this jet is an
    original particle)"""
    parent1::Int

    """Index in _history where second parent of this jet
    was created (NonexistentParent if this jet is an
    original particle); BeamJet if this history entry
    just labels the fact that the jet has recombined
    with the beam)"""
    parent2::Int

    """Index in _history where the current jet is
    recombined with another jet to form its child. It
    is Invalid if this jet does not further
    recombine."""
    child::Int

    """Index in the _jets vector where we will find the
    PseudoJet object corresponding to this jet
    (i.e. the jet created at this entry of the
    history). NB: if this element of the history
    corresponds to a beam recombination, then
    jetp_index=Invalid."""
    jetp_index::Int

    """The distance corresponding to the recombination
       at this stage of the clustering."""
    dij::Float64

    """The largest recombination distance seen
       so far in the clustering history."""
    max_dij_so_far::Float64
end

"""Used for initial particles"""
HistoryElement(jetp_index) = HistoryElement(NonexistentParent, NonexistentParent, Invalid, jetp_index, 0.0, 0.0)


"""
Convienence structure holding all of the relevant parameters for
the jet reconstruction
"""
struct ClusterSequence
    """
    This contains the physical PseudoJets; for each PseudoJet one can find
    the corresponding position in the _history by looking at
    _jets[i].cluster_hist_index()
    """
    jets::Vector{PseudoJet}

    """
    This vector will contain the branching history; for each stage,
    history[i].jetp_index indicates where to look in the _jets
    vector to get the physical PseudoJet.
    """
    history::Vector{HistoryElement}

    """Total energy of the event"""
    Qtot
end

"""Return all inclusive jets of a ClusterSequence with pt > ptmin"""
function inclusive_jets(clusterseq::ClusterSequence, ptmin = 0.0)
    dcut = ptmin * ptmin
    jets_local = Vector{LorentzVectorCyl}(undef, 0)
    # sizehint!(jets_local, length(clusterseq.jets))
    # For inclusive jets with a plugin algorithm, we make no
    # assumptions about anything (relation of dij to momenta,
    # ordering of the dij, etc.)
    # for elt in Iterators.reverse(clusterseq.history)
    for elt in clusterseq.history
        elt.parent2 == BeamJet || continue
        iparent_jet = clusterseq.history[elt.parent1].jetp_index
        jet = clusterseq.jets[iparent_jet]
        if pt2(jet) >= dcut
            push!(jets_local, LorentzVectorCyl(pt(jet), rapidity(jet), phi(jet), mass(jet)))
        end
    end
    jets_local
end

"""Add a new jet's history into the recombination sequence"""
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
