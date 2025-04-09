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
tile_index(tiling_setup, eta::Float64, phi::Float64) = begin
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
tiledjet_set_jetinfo!(jet::TiledJet, clusterseq::ClusterSequence, tiling::Tiling, jets_index, R2, p) = begin
    @inbounds jet.eta = rapidity(clusterseq.jets[jets_index])
    @inbounds jet.phi = phi_02pi(clusterseq.jets[jets_index])
    @inbounds jet.kt2 = pt2(clusterseq.jets[jets_index]) > 1.e-300 ?
                        pt2(clusterseq.jets[jets_index])^p : 1.e300
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
    set_nearest_neighbours!(clusterseq::ClusterSequence, tiling::Tiling, tiledjets::Vector{TiledJet})

This function sets the nearest neighbor information for all jets in the
`tiledjets` vector.

# Arguments
- `clusterseq::ClusterSequence`: The cluster sequence object.
- `tiling::Tiling`: The tiling object.
- `tiledjets::Vector{TiledJet}`: The vector of tiled jets.

# Returns
- `NNs::Vector{TiledJet}`: The vector of nearest neighbor jets.
- `diJ::Vector{Float64}`: The vector of diJ values.

The function iterates over each tile in the `tiling` and sets the nearest
neighbor information for each jet in the tile. It then looks for neighbor jets
in the neighboring tiles and updates the nearest neighbor information
accordingly. Finally, it creates the diJ table and returns the vectors of
nearest neighbor jets and diJ values.

Note: The diJ values are calculated as the kt distance multiplied by R^2.
"""
function set_nearest_neighbours!(clusterseq::ClusterSequence, tiling::Tiling,
                                 tiledjets::Vector{TiledJet})
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

    # Now create the diJ (where J is i's NN) table - remember that
    # we differ from standard normalisation here by a factor of R2
    # (corrected for at the end).
    diJ = similar(clusterseq.jets, Float64)
    NNs = similar(clusterseq.jets, TiledJet)
    for i in eachindex(diJ)
        jetA = tiledjets[i]
        diJ[i] = _tj_diJ(jetA) # kt distance * R^2
        # Our compact diJ table will not be in one-to-one corresp. with non-compact jets,
        # so set up bi-directional correspondence here.
        @inbounds NNs[i] = jetA
        jetA.dij_posn = i
    end
    NNs, diJ
end

"""
    do_ij_recombination_step!(clusterseq::ClusterSequence, jet_i, jet_j, dij, recombine=+)

Perform the bookkeeping associated with the step of recombining jet_i and jet_j
(assuming a distance dij).

# Arguments
- `clusterseq::ClusterSequence`: The cluster sequence object.
- `jet_i`: The index of the first jet to be recombined.
- `jet_j`: The index of the second jet to be recombined.
- `dij`: The distance between the two jets.
- `recombine=+`: The recombination function to be used. Default is addition.

# Returns
- `newjet_k`: The index of the newly created jet.

# Description
This function performs the i-j recombination step in the cluster sequence. It
creates a new jet by recombining the first two jets using the specified
recombination function. The new jet is then added to the cluster sequence. The
function also updates the indices and history information of the new jet and
sorts out the history.
"""
do_ij_recombination_step!(clusterseq::ClusterSequence, jet_i, jet_j, dij, recombine = +) = begin
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

"""
    do_iB_recombination_step!(clusterseq::ClusterSequence, jet_i, diB)

Bookkeeping for recombining a jet with the beam (i.e., finalising the jet) by
adding a step to the history of the cluster sequence.

# Arguments
- `clusterseq::ClusterSequence`: The cluster sequence object.
- `jet_i`: The index of the jet.
- `diB`: The diB value.
"""
do_iB_recombination_step!(clusterseq::ClusterSequence, jet_i, diB) = begin
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
    find_tile_neighbours!(tile_union, jetA, jetB, oldB, tiling)

Find the union of neighbouring tiles of `jetA`, `jetB`, and `oldB` and add them
to the `tile_union`. This established the set of tiles over which searches for
updated and new nearest-neighbours must be run

# Arguments
- `tile_union`: The tile union to which the neighbouring tiles will be added.
- `jetA`: The first jet.
- `jetB`: The second jet.
- `oldB`: The old second jet.
- `tiling`: The tiling information.

# Returns
The number of neighbouring tiles added to the `tile_union`.
"""
function find_tile_neighbours!(tile_union, jetA, jetB, oldB, tiling)
    n_near_tiles = add_untagged_neighbours_to_tile_union(jetA.tile_index,
                                                         tile_union, 0, tiling)
    if isvalid(jetB)
        if jetB.tile_index != jetA.tile_index
            n_near_tiles = add_untagged_neighbours_to_tile_union(jetB.tile_index,
                                                                 tile_union, n_near_tiles,
                                                                 tiling)
        end
        if oldB.tile_index != jetA.tile_index && oldB.tile_index != jetB.tile_index
            n_near_tiles = add_untagged_neighbours_to_tile_union(oldB.tile_index,
                                                                 tile_union, n_near_tiles,
                                                                 tiling)
        end
    end
    n_near_tiles
end

"""
    tiled_jet_reconstruct(particles::AbstractVector{T}; p::Union{Real, Nothing} = -1,
                               algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                               R = 1.0, recombine = +) where {T}

Main jet reconstruction algorithm entry point for reconstructing jets using the
tiled strategy for generic jet type T.

**Note** - if a non-standard recombination is used, it must be defined for
JetReconstruction.PseudoJet, as this struct is used internally.

This code will use the `k_t` algorithm types, operating in `(rapidity, φ)`
space.

It is not necessary to specify both the `algorithm` and the `p` (power) value.
If both are given they must be consistent or an exception is thrown.

## Arguments
- `particles::AbstractVector{T}`: A vector of particles used as input for jet
  reconstruction. T must support methods px, py, pz and energy (defined in the
  JetReconstruction namespace)
- `p::Union{Real, Nothing} = -1`: The power parameter for the jet reconstruction
  algorithm, thus switching between different algorithms.
- `algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing`: The explicit
  jet algorithm to use.
- `R::Float64 = 1.0`: The jet radius parameter for the jet reconstruction
  algorithm.
- `recombine::Function = +`: The recombination function used for combining
  pseudojets.

## Returns
- `Vector{PseudoJet}`: A vector of reconstructed jets.

## Example
```julia
tiled_jet_reconstruct(particles::Vector{LorentzVectorHEP}; p = -1, R = 0.4, recombine = +)
```
"""
function tiled_jet_reconstruct(particles::AbstractVector{T}; p::Union{Real, Nothing} = -1,
                               algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                               R = 1.0, recombine = +) where {T}

    # Check for consistency between algorithm and power
    (p, algorithm) = get_algorithm_power_consistency(p = p, algorithm = algorithm)

    # If we have PseudoJets, we can just call the main algorithm...
    if T == PseudoJet
        # recombination_particles will become part of the cluster sequence, so size it for
        # the starting particles and all N recombinations
        recombination_particles = copy(particles)
        sizehint!(recombination_particles, length(particles) * 2)
    else
        recombination_particles = PseudoJet[]
        sizehint!(recombination_particles, length(particles) * 2)
        for i in eachindex(particles)
            push!(recombination_particles,
                  PseudoJet(px(particles[i]), py(particles[i]), pz(particles[i]),
                            energy(particles[i]), i))
        end
    end

    _tiled_jet_reconstruct(recombination_particles; p = p, R = R, algorithm = algorithm,
                           recombine = recombine)
end

"""
Main jet reconstruction algorithm, using PseudoJet objects
"""

"""
    _tiled_jet_reconstruct(particles::AbstractVector{PseudoJet}; p::Real = -1,
                                algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
                                R = 1.0, recombine = +)

Main jet reconstruction algorithm entry point for reconstructing jets once preprocessing
of data types are done. The algorithm parameter must be consistent with the
power parameter.

## Arguments
- `particles::AbstractVector{PseudoJet}`: A vector of `PseudoJet` particles used as input for jet
  reconstruction.
- `p::Real = -1`: The power parameter for the jet reconstruction algorithm, thus
  switching between different algorithms.
- `R = 1.0`: The jet radius parameter for the jet reconstruction
  algorithm.
- `algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt`: The jet reconstruction
   algorithm to use.
- `recombine::Function = +`: The recombination function used for combining
  pseudojets.

## Returns
- `Vector{PseudoJet}`: A vector of reconstructed jets.

## Example
```julia
tiled_jet_reconstruct(particles::Vector{PseudoJet}; p = 1, R = 1.0, recombine = +)
```
"""
function _tiled_jet_reconstruct(particles::AbstractVector{PseudoJet}; p::Real = -1,
                                algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
                                R = 1.0, recombine = +)
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

    # This will be used quite deep inside loops, so declare it here so that
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

    # ClusterSequence is the struct that holds the state of the reconstruction
    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Tiled, jets, history, Qtot)

    # Tiled jets is a structure that has additional variables for tracking which tile a jet is in
    tiledjets = similar(clusterseq.jets, TiledJet)
    for ijet in eachindex(tiledjets)
        tiledjets[ijet] = TiledJet(ijet)
        tiledjet_set_jetinfo!(tiledjets[ijet], clusterseq, tiling, ijet, R2, p)
    end

    # Now initialise all of the nearest neighbour tiles
    NNs, dij = set_nearest_neighbours!(clusterseq, tiling, tiledjets)

    # Main loop of the reconstruction
    # Each iteration we either merge 2→1 or finalise a jet, so it takes N iterations
    # to complete the reconstruction

    for iteration in 1:N
        # Last slot holds the index of the final valid entry in the
        # compact NNs and dij arrays
        ilast = N - (iteration - 1)
        # Search for the lowest value of min_dij_ijet
        dij_min, ibest = fast_findmin(dij, ilast)
        @inbounds jetA = NNs[ibest]
        jetB = jetA.NN

        # Normalisation
        @fastmath dij_min /= R2

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
            nn = do_ij_recombination_step!(clusterseq, jetA.jets_index, jetB.jets_index,
                                           dij_min, recombine)
            tiledjet_remove_from_tiles!(tiling, jetA)
            oldB = copy(jetB)  # take a copy because we will need it...

            tiledjet_remove_from_tiles!(tiling, jetB)
            tiledjet_set_jetinfo!(jetB, clusterseq, tiling, nn, R2, p) # cause jetB to become _jets[nn]
        #                                  (in addition, registers the jet in the tiling)
        else
            # Jet-beam recombination
            do_iB_recombination_step!(clusterseq, jetA.jets_index, dij_min)
            tiledjet_remove_from_tiles!(tiling, jetA)
            oldB = jetB
        end

        # Find all the neighbour tiles that hold candidate jets for Updates
        n_near_tiles = find_tile_neighbours!(tile_union, jetA, jetB, oldB, tiling)

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
