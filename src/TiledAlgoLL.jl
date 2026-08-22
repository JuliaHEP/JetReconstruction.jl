# Implementation of Tiled Algorithm, using linked list
# This is very similar to FastJet's N2Tiled algorithm
# Original Julia implementation by Philippe Gras,
# ported to this package by Graeme Stewart

using Logging
using Accessors
using LoopVectorization

# Include struct definitions and basic operations
include("TiledAlgoLLStructs.jl")

"""
    _tj_dist(jetA, jetB)

Compute the geometric distance in the (y, ϕ)-plane between two jets in the TiledAlgoLL module.

# Arguments
- `jetA`: The first jet.
- `jetB`: The second jet.

# Returns
The squared distance between `jetA` and `jetB`.

# Examples
"""
_tj_dist(jetA, jetB) = begin
    dphi = π - abs(π - abs(jetA.phi - jetB.phi))
    deta = jetA.eta - jetB.eta
    @muladd dphi * dphi + deta * deta
end

"""
    _tj_diJ(jet)

Compute the dij metric value for a given jet.

# Arguments
- `jet`: The input jet.

# Returns
- The dij value for the jet.

# Example
"""
_tj_diJ(jet) = begin
    kt2 = jet.kt2
    if isvalid(jet.NN) && jet.NN.kt2 < kt2
        kt2 = jet.NN.kt2
    end
    return jet.NN_dist * kt2
end

"""
    tile_index(tiling_setup, eta::Float64, phi::Float64)

Compute the tile index for a given (eta, phi) coordinate.

# Arguments
- `tiling_setup`: The tiling setup object containing the tile size and number of tiles.
- `eta::Float64`: The eta coordinate.
- `phi::Float64`: The phi coordinate.

# Returns
The tile index corresponding to the (eta, phi) coordinate.
"""
function tile_index(tiling_setup, eta::Float64, phi::Float64)
    # Use clamp() to restrict to the correct ranges
    # - eta can be out of range by construction (open ended bins)
    # - phi is protection against bad rounding
    ieta = clamp(1 + unsafe_trunc(Int,
                              (eta - tiling_setup._tiles_eta_min) /
                              tiling_setup._tile_size_eta),
                 1,
                 tiling_setup._n_tiles_eta)
    iphi = clamp(unsafe_trunc(Int, phi / tiling_setup._tile_size_phi), 0,
                 tiling_setup._n_tiles_phi)
    return iphi * tiling_setup._n_tiles_eta + ieta
end

"""
    tiledjet_set_jetinfo!(jet::TiledJet, clusterseq::ClusterSequence, tiling::Tiling, jets_index, R2, p)

Initialise a tiled jet from a PseudoJet (using an index into our ClusterSequence)

Arguments:
- `jet::TiledJet`: The TiledJet object to set the information for.
- `clusterseq::ClusterSequence`: The ClusterSequence object containing the jets.
- `tiling::Tiling`: The Tiling object containing the tile information.
- `jets_index`: The index of the jet in the ClusterSequence.
- `R2`: The jet radius parameter squared.
- `p`: The power to raise the pt2 value to.

This function sets the eta, phi, kt2, jets_index, NN_dist, NN, tile_index, previous, and next fields of the TiledJet object.

Returns:
- `nothing`
"""
function tiledjet_set_jetinfo!(jet::TiledJet,
                               clusterseq::ClusterSequence,
                               tiling::Tiling,
                               jets_index::Int,
                               R2,
                               p)
    real_jet = @inbounds clusterseq.jets[jets_index]

    jet.eta = rapidity(real_jet)
    jet.phi = phi(real_jet)

    jet_pt2 = pt2(real_jet)
    jet.kt2 = jet_pt2 > 1.0e-300 ? jet_pt2^p : 1.0e300

    jet.jets_index = jets_index
    # Initialise NN info as well
    jet.NN_dist = R2
    jet.NN = noTiledJet

    # Find out which tile it belongs to
    jet.tile_index = tile_index(tiling.setup, jet.eta, jet.phi)

    # Insert it into the tile's linked list of jets (at the beginning)
    jet.previous = noTiledJet
    @inbounds jet.next = tiling.tiles[jet.tile_index]

    if isvalid(jet.next)
        jet.next.previous = jet
    end
    @inbounds tiling.tiles[jet.tile_index] = jet
    nothing
end

"""Full scan for nearest neighbours"""

"""
    set_nearest_neighbours!(tiling, tiledjets, NNs, diJ)

Set nearest-neighbour information for all jets in `tiledjets` and fill reusable
nearest-neighbour and distance buffers.

# Arguments
- `tiling::Tiling`: The tiling object.
- `tiledjets::Vector{TiledJet}`: The vector of tiled jets.
- `NNs::Vector{TiledJet}`: Reusable nearest-neighbour storage with the same
  length as `tiledjets`.
- `diJ::Vector{Float64}`: Reusable distance storage with the same length as
  `tiledjets`.

# Returns
- `nothing`

The function iterates over each tile in the `tiling` and sets the nearest
neighbour information for each jet in the tile. It then looks for neighbour jets
in the neighbouring tiles and updates the nearest neighbour information
accordingly. Finally, it fills the supplied nearest-neighbour and diJ tables.

Note: The diJ values are calculated as the kt distance multiplied by R^2.
"""
function set_nearest_neighbours!(tiling::Tiling,
                                 tiledjets::Vector{TiledJet},
                                 NNs::Vector{TiledJet},
                                 diJ::Vector{Float64})
    length(NNs) == length(tiledjets) ||
        throw(ArgumentError("NNs and tiledjets must have the same length"))
    length(diJ) == length(tiledjets) ||
        throw(ArgumentError("diJ and tiledjets must have the same length"))

    # Setup the initial nearest neighbour information
    for tile in tiling.tiles
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
        for rtile_index in rightneighbours(tile.tile_index, tiling)
            for jetA in tile
                for jetB in @inbounds tiling.tiles[rtile_index]
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

    # Now create the diJ (where j is i's NN) table - remember that
    # we differ from standard normalisation here by a factor of R2
    # (corrected for at the end).
    for i in eachindex(diJ)
        @inbounds begin
            jetA = tiledjets[i]

            diJ[i] = _tj_diJ(jetA)
            NNs[i] = jetA
            jetA.dij_posn = i
        end
    end

    return nothing
end

"""
    do_iB_recombination_step!(clusterseq::ClusterSequence, jet_i, diB)

Bookkeeping for recombining a jet with the beam (i.e., finalising the jet) by
adding a step to the history of the cluster sequence.

# Arguments
- `clusterseq::ClusterSequence`: The cluster sequence object.
- `jet_i`: The index of the jet.
- `diB`: The diB value.
"""
function do_iB_recombination_step!(clusterseq::ClusterSequence, jet_i, diB)
    # Recombine the jet with the beam
    add_step_to_history!(clusterseq, clusterseq.jets[jet_i]._cluster_hist_index, BeamJet,
                         Invalid, diB)
end

"""
    add_untagged_neighbours_to_tile_union(center_index, tile_union, n_near_tiles, tiling)

Adds to the vector tile_union the tiles that are in the neighbourhood of the
specified center_index, including itself and whose tagged status are false -
start adding from position n_near_tiles-1, and increase n_near_tiles. When a
neighbour is added its tagged status is set to true.

# Arguments
- `center_index`: The index of the center tile.
- `tile_union`: An array to store the indices of neighbouring tiles.
- `n_near_tiles`: The number of neighbouring tiles.
- `tiling`: The tiling object containing the tile tags.

# Returns
The updated number of near tiles.
"""
function add_untagged_neighbours_to_tile_union(center_index, tile_union, n_near_tiles,
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
    find_tile_neighbours!(tile_union, jetA, jetB, oldB_tile_index::Int, tiling)

Find the union of neighbouring tiles of `jetA`, `jetB`, and `oldB` and add them
to the `tile_union`. This established the set of tiles over which searches for
updated and new nearest-neighbours must be run

# Arguments
- `tile_union`: The tile union to which the neighbouring tiles will be added.
- `jetA`: The first jet.
- `jetB`: The second jet.
- `oldB_tile_index::Int`: The tile index of the old second jet.
- `tiling`: The tiling information.

# Returns
The number of neighbouring tiles added to the `tile_union`.
"""
function find_tile_neighbours!(tile_union, jetA, jetB, oldB_tile_index::Int, tiling)
    n_near_tiles = add_untagged_neighbours_to_tile_union(jetA.tile_index,
                                                         tile_union, 0, tiling)
    if isvalid(jetB)
        if jetB.tile_index != jetA.tile_index
            n_near_tiles = add_untagged_neighbours_to_tile_union(jetB.tile_index,
                                                                 tile_union, n_near_tiles,
                                                                 tiling)
        end
        if oldB_tile_index != jetA.tile_index && oldB_tile_index != jetB.tile_index
            n_near_tiles = add_untagged_neighbours_to_tile_union(oldB_tile_index,
                                                                 tile_union, n_near_tiles,
                                                                 tiling)
        end
    end
    n_near_tiles
end

"""
    tiled_jet_reconstruct(particles::AbstractVector{T};
                          algorithm::JetAlgorithm.Algorithm,
                          p::Union{Real, Nothing} = nothing, R = 1.0,
                          recombine = addjets_escheme, preprocess = preprocess_escheme) where {T}

Main jet reconstruction algorithm entry point for reconstructing jets using the
tiled strategy for generic jet type T.

This code will use the `k_t` algorithm types, operating in `(rapidity, φ)` space.

It is not necessary to specify both the `algorithm` and the `p` (power) value.
If both are given they must be consistent or an exception is thrown.

## Arguments
- `particles::AbstractVector{T}`: A vector of particles used as input for jet
  reconstruction. T must support methods px, py, pz and energy (defined in the
  JetReconstruction namespace).
- `algorithm::JetAlgorithm.Algorithm`: The jet algorithm to use.
- `p::Union{Real, Nothing} = nothing`: The power value used for jet reconstruction.
  Must be specified for GenKt algorithm. Other algorithms will ignore this value.
- `R = 1.0`: The jet radius parameter for the jet reconstruction algorithm.
- `recombine::Function = addjets_escheme`: The recombination function used to combine
  particles into a new jet.
- `preprocess::Function = preprocess_escheme`: A function to preprocess the input particles.

## Returns
- `ClusterSequence`: The resulting `ClusterSequence` object representing the
  reconstructed jets.

## Example
```julia
tiled_jet_reconstruct(particles::Vector{LorentzVectorHEP}; algorithm = JetAlgorithm.GenKt, p = 0.5, R = 0.4)
tiled_jet_reconstruct(particles::Vector{LorentzVectorHEP}; algorithm = JetAlgorithm.AntiKt, R = 0.4)
```
"""
function tiled_jet_reconstruct(particles::AbstractVector{T};
                               algorithm::JetAlgorithm.Algorithm,
                               p::Union{Real, Nothing} = nothing, R = 1.0,
                               recombine = addjets_escheme,
                               preprocess = preprocess_escheme) where {T}

    # Get consistent algorithm power
    p = get_algorithm_power(p = p, algorithm = algorithm)

    recombination_particles = construct_reco_jets(particles,
                                                  PseudoJet,
                                                  preprocess)

    _tiled_jet_reconstruct!(recombination_particles; algorithm = algorithm, p = p, R = R,
                            recombine = recombine)
end

"""
    _tiled_jet_reconstruct!(particles::AbstractVector{PseudoJet};
                           algorithm::JetAlgorithm.Algorithm,
                           p::Real, R = 1.0, recombine = addjets_escheme)

Main jet internal reconstruction algorithm entry point for reconstructing jets
once preprocessing of data types are done. The algorithm parameter must be
consistent with the power parameter.

## Arguments
- `particles::AbstractVector{PseudoJet}`: A vector of `PseudoJet` particles used
  as input for jet reconstruction. This vector must supply the correct
  `cluster_hist_index` values and will be *mutated* as part of the returned
  `ClusterSequence`.
- `algorithm::JetAlgorithm.Algorithm`: The jet reconstruction algorithm to use.
- `p::Real`: The power parameter for the jet reconstruction algorithm, thus
  switching between different algorithms.
- `R = 1.0`: The jet radius parameter for the jet reconstruction algorithm.
- `recombine::Function = addjets_escheme`: The recombination function used for combining
  pseudojets.

## Returns
- `clusterseq`: The resulting `ClusterSequence` object representing the
  reconstructed jets.

## Example
```julia
_tiled_jet_reconstruct!(particles::Vector{PseudoJet}; algorithm = JetAlgorithm.Kt, p = 1, R = 0.4)
```
"""
function _tiled_jet_reconstruct!(particles::AbstractVector{PseudoJet};
                                 algorithm::JetAlgorithm.Algorithm,
                                 p::Real,
                                 R = 1.0,
                                 recombine = addjets_escheme,
                                 scratch::Union{Nothing, TiledScratch} = nothing,
                                 history_buffer::Union{Nothing, Vector{HistoryElement}} = nothing)
    # Bounds
    N::Int = length(particles)

    # Extremely odd - having these @debug statements present causes a performance
    # degradation of ~20μs per event on my M2 mac (12%!), even when no debugging is used
    # so they need to be completely commented out...
    #
    # There are a few reports of this in, e.g., https://github.com/JuliaLang/julia/issues/28147
    # It does seem to have improved, but it's far from perfect!
    # @debug "Initial particles: $(N)"

    # Algorithm parameters
    R2::Float64 = R * R
    p = (round(p) == p) ? Int(p) : p # integer p if possible

    # Normalise the optional storage once. The owning path uses temporary
    # buffers, while the workspace path supplies reusable ones; reconstruction
    # below is identical in both cases.
    scratch = isnothing(scratch) ? TiledScratch() : scratch
    ensure_length!(scratch, N)

    tile_union = scratch.tile_union
    eta_local = scratch.eta
    tiledjets = scratch.tiledjets

    if isnothing(history_buffer)
        # Preserve the owning path's eager complete-history allocation.
        history_buffer = Vector{HistoryElement}(undef, N)
        sizehint!(history_buffer, 2 * N)
    end
    history, Qtot = initial_history!(history_buffer, particles)

    # Now get the tiling setup
    for ijet in 1:N
        @inbounds eta_local[ijet] = rapidity(particles[ijet])
    end

    setup = setup_tiling(eta_local, R)

    tiling = get_tiling!(scratch, setup)

    # ClusterSequence is the struct that holds the state of the reconstruction
    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Tiled, particles, history,
                                 Qtot)

    # Tiled jets is a structure that has additional variables for tracking which tile a jet is in
    for ijet in 1:N
        @inbounds begin
            jet = tiledjets[ijet]
            # A previous event may have repurposed this node during merging.
            # Restore its current event-local initial identity.
            jet.id = ijet

            tiledjet_set_jetinfo!(jet, clusterseq, tiling, ijet, R2, p)
        end
    end

    # Now initialise all of the nearest neighbour tiles
    NNs = scratch.NNs
    dij = scratch.dij
    set_nearest_neighbours!(tiling, tiledjets, NNs, dij)

    # Main loop of the reconstruction
    # Each iteration we either merge 2→1 or finalise a jet, so it takes N iterations
    # to complete the reconstruction

    for iteration in 1:N
        # Last slot holds the index of the final valid entry in the
        # compact NNs and diJ arrays
        ilast = N - (iteration - 1)
        # Search for the lowest value of min_dij_ijet
        dij_min, ibest = fast_findmin(dij, ilast)
        @inbounds jetA = NNs[ibest]
        jetB = jetA.NN

        # Normalisation
        @fastmath dij_min /= R2

        # @debug "Iteration $(iteration): dij_min $(dij_min); jetA $(jetA.id), jetB $(jetB.id)"
        oldB_tile_index = 0
        if isvalid(jetB)
            # Jet-jet recombination
            # If necessary relabel A & B to ensure jetB < jetA, that way if
            # the larger of them == newtail then that ends up being jetA and
            # the new jet that is added as jetB is inserted in a position that
            # has a future!
            if jetA.id < jetB.id
                jetA, jetB = jetB, jetA
            end

            # Only the old tile index is needed below. Copying this mutable
            # TiledJet would allocate once for every jet-jet recombination.
            oldB_tile_index = jetB.tile_index

            # Recombine jetA and jetB into the new jet
            real_jetA = clusterseq.jets[jetA.jets_index]
            real_jetB = clusterseq.jets[jetB.jets_index]
            newjet = recombine(real_jetA, real_jetB;
                               cluster_hist_index = length(clusterseq.history) + 1)
            push!(clusterseq.jets, newjet)
            newjet_k = length(clusterseq.jets)
            add_step_to_history!(clusterseq,
                                 minmax(real_jetA._cluster_hist_index,
                                        real_jetB._cluster_hist_index)...,
                                 newjet_k, dij_min)

            tiledjet_remove_from_tiles!(tiling, jetA)

            tiledjet_remove_from_tiles!(tiling, jetB)
            # Move jetB to be jets[newjet_k] and register the new jet in the tiling
            tiledjet_set_jetinfo!(jetB, clusterseq, tiling, newjet_k, R2, p)
        else
            # Jet-beam recombination
            do_iB_recombination_step!(clusterseq, jetA.jets_index, dij_min)
            tiledjet_remove_from_tiles!(tiling, jetA)
        end

        # Find all the neighbour tiles that hold candidate jets for Updates
        n_near_tiles = find_tile_neighbours!(tile_union, jetA, jetB, oldB_tile_index,
                                             tiling)

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
                    jetI.NN = noTiledJet

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
                        jetB.NN = jetI
                    end
                end # isvalid(jetB)
            end #next jetI
        end #next itile

        # finally, register the updated kt distance for B
        if isvalid(jetB)
            @inbounds dij[jetB.dij_posn] = _tj_diJ(jetB)
        end
    end
    clusterseq
end

"""
Run full-semantics N2Tiled reconstruction using workspace-owned storage.

The caller must ensure exclusive ownership of the workspace.
"""
function _n2tiled_reconstruct_with_workspace!(workspace::N2TiledWorkspace,
                                              particles::AbstractVector;
                                              algorithm::JetAlgorithm.Algorithm,
                                              p::Union{Real, Nothing} = nothing,
                                              R = 1.0,
                                              recombine = addjets_escheme,
                                              preprocess = preprocess_escheme)
    resolved_power = get_algorithm_power(p = p,
                                         algorithm = algorithm)

    jets = construct_reco_jets!(workspace.jets,
                                particles,
                                preprocess)

    return _tiled_jet_reconstruct!(jets;
                                   algorithm = algorithm,
                                   p = resolved_power,
                                   R = R,
                                   recombine = recombine,
                                   scratch = workspace.scratch,
                                   history_buffer = workspace.history)
end

"""
    with_n2tiled_reconstruction(
        f,
        workspace,
        particles;
        algorithm,
        p=nothing,
        R=1.0,
        recombine=addjets_escheme,
        preprocess=preprocess_escheme,
    )

Run full N2Tiled reconstruction using workspace-owned reusable storage, then
invoke `f` with the resulting borrowed ClusterSequence.

# Borrowed lifetime

Inside the callback:

- `clusterseq.jets === workspace.jets`;
- `clusterseq.history === workspace.history`;
- complete clustering history is available;
- constituent and exclusive-jet queries remain available.

The ClusterSequence, its jets vector, and its history vector must not be retained
for use after this workspace processes another event.

Copy any data that needs to outlive the callback.

The same workspace must not be used concurrently or reentrantly.
Use ordinary `jet_reconstruct` or `tiled_jet_reconstruct` when an independently
owned ClusterSequence is required.
"""
function with_n2tiled_reconstruction(f::F,
                                     workspace::N2TiledWorkspace,
                                     particles::AbstractVector{T};
                                     algorithm::JetAlgorithm.Algorithm,
                                     p::Union{Real, Nothing} = nothing,
                                     R = 1.0,
                                     recombine = addjets_escheme,
                                     preprocess = preprocess_escheme) where {F, T}
    clusterseq = _n2tiled_reconstruct_with_workspace!(workspace,
                                                      particles;
                                                      algorithm = algorithm,
                                                      p = p,
                                                      R = R,
                                                      recombine = recombine,
                                                      preprocess = preprocess)

    return f(clusterseq)
end
