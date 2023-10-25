# Implementation of Tiled Algorithm, using linked list
# This is very similar to FastJet's N2Tiled algorithm
# Original Julia implementation by Philippe Gras,
# ported to this package by Graeme Stewart

using Logging
using Accessors
using LoopVectorization

# Include struct definitions and basic operations
include("TiledAlglLLStructs.jl")

"""
Initialise the clustering history in a standard way,
Takes as input the list of stable particles as input
Returns the history and the total event energy.
"""
function initial_history(particles)
    # reserve sufficient space for everything
    history = Vector{HistoryElement}(undef, length(particles))
    sizehint!(history, 2 * length(particles))

    Qtot::Float64 = 0

    for i in eachindex(particles)
        history[i] = HistoryElement(i)

        # get cross-referencing right from PseudoJets
        particles[i]._cluster_hist_index = i

        # determine the total energy in the event
        Qtot += particles[i].E
    end
    history, Qtot
end

"""Computes distance in the (eta,phi)-plane
between two jets."""
_tj_dist(jetA, jetB) = begin
    dphi = π - abs(π - abs(jetA.phi - jetB.phi))
    deta = jetA.eta - jetB.eta
    return dphi * dphi + deta * deta
end

_tj_diJ(jet) = begin
    kt2 = jet.kt2
    if isvalid(jet.NN) && jet.NN.kt2 < kt2
        kt2 = jet.NN.kt2
    end
    return jet.NN_dist * kt2
end


"""Return the tile index corresponding to the given eta,phi point"""
tile_index(tiling_setup, eta::Float64, phi::Float64) = begin
    # Use clamp() to restrict to the correct ranges
    # - eta can be out of range by construction (open end bins)
    # - phi is protection against bad rounding
    ieta = clamp(1 + unsafe_trunc(Int, (eta - tiling_setup._tiles_eta_min) / tiling_setup._tile_size_eta), 1, tiling_setup._n_tiles_eta)
    iphi = clamp(unsafe_trunc(Int, phi / tiling_setup._tile_size_phi), 0, tiling_setup._n_tiles_phi)
    return iphi * tiling_setup._n_tiles_eta + ieta
end


"""Initialise a tiled jet from a PseudoJet (using an index into our ClusterSequence)"""
tiledjet_set_jetinfo!(jet::TiledJet, clusterseq::ClusterSequence, jets_index, R2) = begin
    @inbounds jet.eta = rapidity(clusterseq.jets[jets_index])
    @inbounds jet.phi = phi_02pi(clusterseq.jets[jets_index])
    @inbounds jet.kt2 = pt2(clusterseq.jets[jets_index]) > 1.e-300 ? 1.0 / pt2(clusterseq.jets[jets_index]) : 1.e300
    jet.jets_index    = jets_index
    # Initialise NN info as well
    jet.NN_dist = R2
    jet.NN      = noTiledJet

    # Find out which tile it belonds to
    jet.tile_index = tile_index(clusterseq.tiling.setup, jet.eta, jet.phi)

    # Insert it into the tile's linked list of jets (at the beginning)
    jet.previous = noTiledJet
    @inbounds jet.next = clusterseq.tiling.tiles[jet.tile_index]
    if isvalid(jet.next)
        jet.next.previous = jet
    end
    @inbounds clusterseq.tiling.tiles[jet.tile_index] = jet
    nothing
end


"""Full scan for nearest neighbours"""
function set_nearest_neighbours!(clusterseq::ClusterSequence, tiledjets::Vector{TiledJet})
    # Setup the initial nearest neighbour information
    for tile in clusterseq.tiling.tiles
        isvalid(tile) || continue
        for jetA in tile
            for jetB in tile
                if jetB == jetA
                    break
                end
                dist = _tj_dist(jetA, jetB)
                if (dist < jetA.NN_dist)
                    jetA.NN_dist = dist
                    jetA.NN = jetB
                end
                if dist < jetB.NN_dist
                    jetB.NN_dist = dist
                    jetB.NN = jetA
                end
            end
        end

        # Look for neighbour jets n the neighbour tiles
        for rtile_index in rightneighbours(tile.tile_index, clusterseq.tiling)
            for jetA in tile
                for jetB in @inbounds clusterseq.tiling.tiles[rtile_index]
                    dist = _tj_dist(jetA, jetB)
                    if (dist < jetA.NN_dist)
                        jetA.NN_dist = dist
                        jetA.NN = jetB
                    end
                    if dist < jetB.NN_dist
                        jetB.NN_dist = dist
                        jetB.NN = jetA
                    end
                end
            end
            # No need to do it for LH tiles, since they are implicitly done
            # when we set NN for both jetA and jetB on the RH tiles.
        end
    end

    # Now create the diJ (where J is i's NN) table - remember that
    # we differ from standard normalisation here by a factor of R2
    # (corrected for at the end).
    diJ = similar(clusterseq.jets, Float64)
    NNs = similar(clusterseq.jets, TiledJet)
    for i in eachindex(diJ)
        jetA = tiledjets[i]
        diJ[i] = _tj_diJ(jetA) # kt distance * R^2
        # our compact diJ table will not be in one-to-one corresp. with non-compact jets,
        # so set up bi-directional correspondence here.
        @inbounds NNs[i] = jetA
        jetA.dij_posn = i
    end
    NNs, diJ
end

"""Carries out the bookkeeping associated with the step of recombining
jet_i and jet_j (assuming a distance dij) and returns the index
of the recombined jet, newjet_k."""
do_ij_recombination_step!(clusterseq::ClusterSequence, jet_i, jet_j, dij, recombine=+) = begin
    # Create the new jet by recombining the first two with
    # the E-scheme
    push!(clusterseq.jets, recombine(clusterseq.jets[jet_i], clusterseq.jets[jet_j]))

    # Get its index and the history index
    newjet_k = length(clusterseq.jets)
    newstep_k = length(clusterseq.history) + 1

    # And provide jet with this info
    clusterseq.jets[newjet_k]._cluster_hist_index = newstep_k

    # Finally sort out the history
    hist_i = clusterseq.jets[jet_i]._cluster_hist_index
    hist_j = clusterseq.jets[jet_j]._cluster_hist_index

    add_step_to_history!(clusterseq, minmax(hist_i, hist_j)...,
        newjet_k, dij)

    newjet_k
end

"""Carries out the bookkeeping associated with the step of recombining
jet_i with the beam (i.e., finalising the jet)"""
do_iB_recombination_step!(clusterseq::ClusterSequence, jet_i, diB) = begin
    # Recombine the jet with the beam
    add_step_to_history!(clusterseq, clusterseq.jets[jet_i]._cluster_hist_index, BeamJet,
        Invalid, diB)
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
        throw(ErrorException("Internal error. Trying to recombine an object that has previsously been recombined."))
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

"""Adds to the vector tile_union the tiles that are in the neighbourhood
of the specified tile_index, including itself and whose tagged status are
false ---start adding from position n_near_tiles-1, and increase n_near_tiles as
you go along. When a neighbour is added its tagged status is set to true.

Returns the updated number of near_tiles."""
function add_untagged_neighbours_to_tile_union(center_index,
    tile_union, n_near_tiles,
    tiling)
    for tile_index in surrounding(center_index, tiling)
        @inbounds if !tiling.tags[tile_index]
            n_near_tiles += 1
            tile_union[n_near_tiles] = tile_index
            tiling.tags[tile_index] = true
        else
        end
    end
    n_near_tiles
end

"""
Establish the set of tiles over which we are going to
have to run searches for updated and new nearest-neighbours --
basically a combination of vicinity of the tiles of the two old
jets and the new jet.

Updates tile_union and returns n_near_tiles
"""
function find_tile_neighbours!(tile_union, jetA, jetB, oldB, tiling)
    n_near_tiles = add_untagged_neighbours_to_tile_union(jetA.tile_index,
        tile_union, 0, tiling)
    if isvalid(jetB)
        if jetB.tile_index != jetA.tile_index
            n_near_tiles = add_untagged_neighbours_to_tile_union(jetB.tile_index,
                tile_union, n_near_tiles, tiling)
        end
        if oldB.tile_index != jetA.tile_index && oldB.tile_index != jetB.tile_index
            n_near_tiles = add_untagged_neighbours_to_tile_union(oldB.tile_index,
                tile_union, n_near_tiles, tiling)
        end
    end
    n_near_tiles
end

"""Find the lowest value in the array, returning the value and the index"""
find_lowest(dij, n) = begin
    best = 1
    @inbounds dij_min = dij[1]
    @turbo for here in 2:n
        newmin = dij[here] < dij_min
        best = newmin ? here : best
        dij_min = newmin ? dij[here] : dij_min
    end
    # @assert argmin(dij[1:n]) == best
    dij_min, best
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
            # push!(jets_local, jet)
        end
    end
    jets_local
end


"""
Main jet reconstruction algorithm entry point for generic data types

`particles` must support methods px, py, pz and energy (N.B. these must be in the
JetReconstruction namespace). In particular Cartesian LorentzVector structs can
be used.

If a non-standard recombination is used, it must be defined for
JetReconstruction.PseudoJet, as this struct is used internally.
"""
function tiled_jet_reconstruct(particles::Vector{T}; p = -1, R = 1.0, recombine = +, ptmin = 0.0) where {T}
    # Here we need to populate the vector of PseudoJets that are the internal
    # EDM for the main algorithm, then we call the reconstruction
    pseudojets = Vector{PseudoJet}(undef, length(particles))
    for (i, particle) in enumerate(particles)
        pseudojets[i] = PseudoJet(px(particle), py(particle),
            pz(particle), energy(particle))
    end
    tiled_jet_reconstruct(pseudojets, p = p, R = R, recombine = recombine, ptmin = ptmin)
end

"""
Main jet reconstruction algorithm, using PseudoJet objects
"""
function tiled_jet_reconstruct(particles::Vector{PseudoJet}; p = -1, R = 1.0, recombine = +, ptmin = 0.0)
    # Bounds
    N::Int = length(particles)
    # @debug "Initial particles: $(N)"

    # Algorithm parameters
    R2::Float64 = R * R
    p = (round(p) == p) ? Int(p) : p # integer p if possible

    # This will be used quite deep inside loops, but declare it here so that
    # memory (de)allocation gets done only once
    tile_union = Vector{Int}(undef, 3 * _n_tile_neighbours)

    # Container for pseudojets, sized for all initial particles, plus all of the
    # merged jets that can be created during reconstruction
    jets = PseudoJet[]
    sizehint!(jets, N * 2)
    resize!(jets, N)

    # Copy input data into the jets container
    copyto!(jets, particles)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(jets)

    # Now get the tiling setup
    _eta = Vector{Float64}(undef, length(particles))
    for ijet in 1:length(particles)
        _eta[ijet] = rapidity(particles[ijet])
    end

    tiling = Tiling(setup_tiling(_eta, R))

    # ClusterSequence is a convenience struct that holds the state of the reconstruction
    clusterseq = ClusterSequence(jets, history, tiling, Qtot)

    # Tiled jets is a structure that has additional variables for tracking which tile a jet is in
    tiledjets = similar(clusterseq.jets, TiledJet)
    for ijet in eachindex(tiledjets)
        tiledjets[ijet] = TiledJet(ijet)
        tiledjet_set_jetinfo!(tiledjets[ijet], clusterseq, ijet, R2)
    end

    # Now initalise all of the nearest neighbour tiles
    NNs, dij = set_nearest_neighbours!(clusterseq, tiledjets)

    # Main loop of the reconstruction
    # Each iteration we either merge 2→1 or finalise a jet, so it takes N iterations
    # to complete the reconstruction

    for iteration in 1:N
        # Last slot holds the index of the final valid entry in the
        # compact NNs and dij arrays
        ilast = N - (iteration - 1)
        # Search for the lowest value of min_dij_ijet
        dij_min, ibest = find_lowest(dij, ilast)
        @inbounds jetA = NNs[ibest]
        jetB = jetA.NN

        # Normalisation
        dij_min *= R2

        # @debug "Iteration $(iteration): dij_min $(dij_min); jetA $(jetA.id), jetB $(jetB.id)"

        if isvalid(jetB)
            # Jet-jet recombination
            # If necessary relabel A & B to ensure jetB < jetA, that way if
            # the larger of them == newtail then that ends up being jetA and
            # the new jet that is added as jetB is inserted in a position that
            # has a future!
            if jetA.id < jetB.id
                jetA, jetB = jetB, jetA
            end

            # Recombine jetA and jetB and retrieves the new index, nn
            nn = do_ij_recombination_step!(clusterseq, jetA.jets_index, jetB.jets_index, dij_min, recombine)
            tiledjet_remove_from_tiles!(clusterseq.tiling, jetA)
            oldB = copy(jetB)  # take a copy because we will need it...

            tiledjet_remove_from_tiles!(clusterseq.tiling, jetB)
            tiledjet_set_jetinfo!(jetB, clusterseq, nn, R2) # cause jetB to become _jets[nn]
        #                                  (in addition, registers the jet in the tiling)
        else
            # Jet-beam recombination
            do_iB_recombination_step!(clusterseq, jetA.jets_index, dij_min)
            tiledjet_remove_from_tiles!(clusterseq.tiling, jetA)
            oldB = jetB
        end

        # Find all the neighbour tiles that hold candidate jets for Updates
        n_near_tiles = find_tile_neighbours!(tile_union, jetA, jetB, oldB, clusterseq.tiling)

        # Firstly compactify the diJ by taking the last of the diJ and copying
        # it to the position occupied by the diJ for jetA
        @inbounds NNs[ilast].dij_posn = jetA.dij_posn
        @inbounds dij[jetA.dij_posn] = dij[ilast]
        @inbounds NNs[jetA.dij_posn] = NNs[ilast]

        # Initialise jetB's NN distance as well as updating it for
        # other particles.
        # Run over all tiles in our union
        for itile in 1:n_near_tiles
            @inbounds tile = tiling.tiles[@inbounds tile_union[itile]] #TAKES 5μs
            @inbounds tiling.tags[tile_union[itile]] = false # reset tag, since we're done with unions

            isvalid(tile) || continue #Probably not required

            # run over all jets in the current tile
            for jetI in tile
                # see if jetI had jetA or jetB as a NN -- if so recalculate the NN
                if jetI.NN == jetA || (jetI.NN == jetB && isvalid(jetB))
                    jetI.NN_dist = R2
                    jetI.NN      = noTiledJet

                    # now go over tiles that are neighbours of I (include own tile)
                    for near_tile_index in surrounding(tile.tile_index, tiling)
                        # and then over the contents of that tile
                        for jetJ in @inbounds tiling.tiles[near_tile_index]
                            # Sometimes jetJ comes out as invalid...?
                            dist = _tj_dist(jetI, jetJ)
                            if dist < jetI.NN_dist && jetJ != jetI
                                jetI.NN_dist = dist
                                jetI.NN = jetJ
                            end
                        end # next jetJ
                    end # next near_tile
                    dij[jetI.dij_posn] = _tj_diJ(jetI) # update diJ kt-dist
                end #jetI.NN == jetA || (jetI.NN == jetB && !isnothing(jetB))

                # check whether new jetB is closer than jetI's current NN and
                # if jetI is closer than jetB's current (evolving) nearest
                # neighbour. Where relevant update things.
                if isvalid(jetB)
                    dist = _tj_dist(jetI, jetB)
                    if dist < jetI.NN_dist
                        if jetI != jetB
                            jetI.NN_dist = dist
                            jetI.NN = jetB
                            dij[jetI.dij_posn] = _tj_diJ(jetI) # update diJ...
                        end
                    end
                    if dist < jetB.NN_dist && jetI != jetB
                        jetB.NN_dist = dist
                        jetB.NN      = jetI
                    end
                end # isvalid(jetB)
            end #next jetI
        end #next itile

        # finally, register the updated kt distance for B
        if isvalid(jetB)
            @inbounds dij[jetB.dij_posn] = _tj_diJ(jetB)
        end
    end
    inclusive_jets(clusterseq, ptmin), clusterseq.history
end
