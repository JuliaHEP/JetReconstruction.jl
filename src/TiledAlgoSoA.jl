### Tiled Jet Reconstruction, SoA

using Logging

"""
Do a complete scan over tiles to find all nearest neighbour distances
"""
function find_all_nearest_neighbours!(tile_jets::Array{TiledJetSoA, 2}, tiling_setup::TilingDef,
	flat_jets::FlatJetSoA, R2::AbstractFloat)
	# We march over all tiles, evaluating the distance to
	# - each jet in this tile
	# - each jet in the rightmost neighbour tiles
	# As we compare jet-to-jet, in both directions, the rightmost tiles ensure that the
	# whole space is swept
	@inbounds for itile in eachindex(tile_jets)
		tile = tile_jets[itile]
		if (tile._size == 0)
			continue
		end
		@inbounds for ijet in 1:tile._size
			# Could we do this in a broadcast way...? Would it be faster?
			@inbounds for jjet in ijet+1:tile._size
				_dist = geometric_distance(tile._eta[ijet], tile._phi[ijet],
					tile._eta[jjet], tile._phi[jjet])
				if (_dist < tile._nndist[ijet])
					tile._nndist[ijet] = _dist
					set_nn!(tile._nn[ijet], itile, jjet)
				end
				if (_dist < tile._nndist[jjet])
					tile._nndist[jjet] = _dist
					set_nn!(tile._nn[jjet], itile, ijet)
				end
			end
			@inbounds for jtile in tile._righttiles
				tile2 = tile_jets[jtile]
				@inbounds for jjet in 1:tile2._size
					_dist = geometric_distance(tile._eta[ijet], tile._phi[ijet],
						tile2._eta[jjet], tile2._phi[jjet])
					if (_dist < tile._nndist[ijet])
						tile._nndist[ijet] = _dist
						set_nn!(tile._nn[ijet], jtile, jjet)
					end
					if (_dist < tile2._nndist[jjet])
						tile2._nndist[jjet] = _dist
						set_nn!(tile2._nn[jjet], itile, ijet)
					end
				end
			end
		end
	end

	# Now calculate the dij distances
	min_dij = 1e20
	min_dij_itile = 0
	min_dij_ijet = 0
	@inbounds for (itile, tile) in enumerate(tile_jets)
		@inbounds for ijet in 1:tile._size
			if valid_nn(tile._nn[ijet])
				tile._dij[ijet] = tile._nndist[ijet] *
								  min(tile._kt2[ijet], tile_jets[tile._nn[ijet]._itile]._kt2[tile._nn[ijet]._ijet])
			else
				tile._dij[ijet] = tile._nndist[ijet] * tile._kt2[ijet]
			end
			# And as this is a complete scan, find the minimum dij as well
			if tile._dij[ijet] < min_dij
				min_dij = tile._dij[ijet]
				min_dij_itile = itile
				min_dij_ijet = ijet
			end
		end
	end
	min_dij_itile, min_dij_ijet, min_dij
end

"""
Scan a tile and all of its neighbours for jet distances
"""
function scan_neighbors!(tile_jets::Array{TiledJetSoA, 2}, jet_tile_index::TiledSoACoord, R2::AbstractFloat)
    itile = jet_tile_index._itile
	tile = tile_jets[itile]
	ijet = jet_tile_index._ijet
    # println("$(tile) $(ijet) $(tile._size)")
	@inbounds for jjet in 1:tile._size
		if jjet == ijet
			continue
		end
		_dist = geometric_distance(tile._eta[ijet], tile._phi[ijet],
			tile._eta[jjet], tile._phi[jjet])
		if (_dist < tile._nndist[ijet])
			tile._nndist[ijet] = _dist
			set_nn!(tile._nn[ijet], itile, jjet)
		end
		if (_dist < tile._nndist[jjet])
			tile._nndist[jjet] = _dist
			set_nn!(tile._nn[jjet], itile, ijet)
			tile._dij[jjet] = tile._nndist[jjet] * min(tile._kt2[ijet], tile._kt2[jjet])
		end
	end
	@inbounds for jtile in tile._nntiles
        tile_j = tile_jets[jtile]
		@inbounds for jjet in 1:tile_j._size
			_dist = geometric_distance(tile._eta[ijet], tile._phi[ijet],
				tile_j._eta[jjet], tile_j._phi[jjet])
			if (_dist < tile._nndist[ijet])
				tile._nndist[ijet] = _dist
				set_nn!(tile._nn[ijet], jtile, jjet)
			end
			if (_dist < tile_j._nndist[jjet])
				tile_j._nndist[jjet] = _dist
				set_nn!(tile_j._nn[jjet], itile, ijet)
				tile_j._dij[jjet] = tile_j._nndist[jjet] * min(tile._kt2[ijet], tile_j._kt2[jjet])
			end
		end
	end
	if tile._nndist[ijet] != R2
		tile._dij[ijet] = tile._nndist[ijet] * min(tile._kt2[ijet],
			tile_jets[tile._nn[ijet]._itile]._kt2[tile._nn[ijet]._ijet])
	end

end

"""Dump NN and dij distances for the current tiled jets, for debugging"""
function get_nn_str(tile_jets::Array{TiledJetSoA, 2})
	jet_nn_strings = Vector{Tuple{Int64, String}}(undef, 0)
	for itile in eachindex(tile_jets)
		tile = tile_jets[itile]
		for ijet in 1:tile._size
			push!(jet_nn_strings, (tile._index[ijet], "Jet $(tile._index[ijet]), NN is $(nnindex(tile_jets, itile, ijet)), geo $(tile._nndist[ijet]) dij $(tile._dij[ijet])\n"))
		end
	end
	sort!(jet_nn_strings)
	nn_str = ""
	for (ijet, jet_str) in jet_nn_strings
		nn_str *= jet_str
	end
	nn_str
end

"""Validate distance calculations, for debugging"""
function validate_distances(tile_jets::Array{TiledJetSoA, 2}, flat_jets, ijetA::Int, ijetB::Int)
	itileA = itjetA = itileB = itjetB = 0
	for itile in eachindex(tile_jets)
		tile = tile_jets[itile]
		for ijet in 1:tile._size
			if tile._index[ijet] == ijetA
				itileA = itile
				itjetA = ijet
			elseif tile._index[ijet] == ijetB
				itileB = itile
				itjetB = ijet
			end
		end
	end
	dist_string = "$(ijetA) $(flat_jets._eta[ijetA]) $(flat_jets._phi[ijetA]) $(flat_jets._kt2[ijetA]) ($(itileA),$(itjetA))\n"
	dist_string *= "$(ijetB) $(flat_jets._eta[ijetB]) $(flat_jets._phi[ijetB]) $(flat_jets._kt2[ijetB]) ($(itileB),$(itjetB))\n"
	dist_string *= "$(tile_jets[itileA]._eta[itjetA]) $(tile_jets[itileA]._phi[itjetA])\n"
	dist_string *= "$(tile_jets[itileB]._eta[itjetB]) $(tile_jets[itileB]._phi[itjetB])\n"
	dist = geometric_distance(tile_jets[itileA]._eta[itjetA], tile_jets[itileA]._phi[itjetA], tile_jets[itileB]._eta[itjetB], tile_jets[itileB]._phi[itjetB])
	dist_string *= "Geo $(dist)"
	dist_string
end

"""
Tiled jet reconstruction
"""
function tiled_jet_reconstruct_soa(objects::AbstractArray{T}; p = -1, R = 1.0, recombine = +) where T
	# bounds
	N::Int = length(objects)
	@debug "Initial particles: $(N)"

	# params
	_R2::Float64 = R * R
	_p = (round(p) == p) ? Int(p) : p # integer p if possible
	ap = abs(_p) # absolute p

	# Input data
	_objects = copy(objects)
	sizehint!(_objects, N * 2)
	_kt2 = (JetReconstruction.pt.(_objects) .^ 2) .^ _p
	sizehint!(_kt2, N * 2)
	_phi = JetReconstruction.phi.(_objects)
	sizehint!(_phi, N * 2)
	_eta = JetReconstruction.eta.(_objects)
	sizehint!(_eta, N * 2)
	_index = collect(1:N) # Initial jets are just numbered 1:N
	sizehint!(_index, N * 2)

	# returned values
	jets = T[] # result
	_sequences = Vector{Int}[[x] for x in 1:N]

	flat_jets = FlatJetSoA(N, _kt2, _eta, _phi, _index)

	# Tiling
	tiling_setup = setup_tiling(_eta, R)
	tile_jets = Array{TiledJetSoA, 2}(undef, tiling_setup._n_tiles_eta, tiling_setup._n_tiles_phi)

	# Populate tiles, from the initial particles
	populate_tiles!(tile_jets, tiling_setup, flat_jets, _R2)

	# Setup initial nn, nndist and dij values
	min_dij_itile, min_dij_ijet, min_dij = find_all_nearest_neighbours!(tile_jets, tiling_setup, flat_jets, _R2)

	# Move some variables outside the loop, to avoid per-loop allocations
	itouched_tiles = Set{Int}()
	sizehint!(itouched_tiles, 12)
	tainted_slots = Set{TiledSoACoord}()
	sizehint!(tainted_slots, 4)

	# At each iteration we either merge two jets to one, or finalise a jet
	# Thus each time we lose one jet, and it therefore takes N iterations to complete
	# the algorithm
	for iteration in 1:N
        @debug "Iteration $(iteration)"

		# For the first iteration the nearest neighbour is known
		if iteration > 1
			# Now find again the new nearest dij jets
			min_dij = 1.0e20
			min_dij_itile = 0
			min_dij_ijet = 0
			@inbounds for itile in eachindex(tile_jets)
				@inbounds for ijet in 1:tile_jets[itile]._size
					if tile_jets[itile]._dij[ijet] < min_dij
						min_dij_itile = itile
						min_dij_ijet = ijet
						min_dij = tile_jets[itile]._dij[ijet]
					end
				end
			end
		end

        @debug "$(min_dij) at ($(min_dij_itile), $(min_dij_ijet)) $(tile_jets[min_dij_itile]._index[min_dij_ijet]) -> $(tile_jets[min_dij_itile]._nn[min_dij_ijet])"
		# Is this a merger or a final jet?
		if tile_jets[min_dij_itile]._nn[min_dij_ijet]._itile == 0
			# Final jet
            jet_merger = false
			index_tile_jetA = TiledSoACoord(min_dij_itile, min_dij_ijet)
			index_jetA = tile_jets[min_dij_itile]._index[min_dij_ijet]
			empty!(tainted_slots)
			push!(tainted_slots, index_tile_jetA)
            push!(jets, _objects[index_jetA])
			push!(_sequences[index_jetA], 0)
			@debug "Finalise jet $(tile_jets[min_dij_itile]._index[min_dij_ijet]) ($(_sequences[index_jetA])) $(JetReconstruction.pt(_objects[index_jetA]))"
			push!(tainted_slots, remove_jet!(tile_jets, index_tile_jetA))
		else
			# Merge two jets
            jet_merger = true
			index_tile_jetA = TiledSoACoord(min_dij_itile, min_dij_ijet)
			index_tile_jetB = tile_jets[min_dij_itile]._nn[min_dij_ijet]
			index_jetA = tile_jets[min_dij_itile]._index[min_dij_ijet]
			index_jetB = nnindex(tile_jets, min_dij_itile, min_dij_ijet)
			@debug "Merge jets $(index_jetA) ($(index_tile_jetA)) and $(index_jetB) ($(index_tile_jetB))"
			merged_jet = recombine(_objects[index_jetA], _objects[index_jetB])

            # If A and B are in the same tile, ensure that A is the earlier slot
            # so that slots are filled up correctly
            if (index_tile_jetA._itile == index_tile_jetB._itile) && (index_tile_jetA._ijet > index_tile_jetB._ijet)
                index_tile_jetA, index_tile_jetB = index_tile_jetB, index_tile_jetA
                index_jetA, index_jetB = index_jetB, index_jetA
            end

			push!(_objects, merged_jet)
			push!(flat_jets._index, length(_objects))
			push!(flat_jets._phi, JetReconstruction.phi(merged_jet))
			push!(flat_jets._eta, JetReconstruction.eta(merged_jet))
			push!(flat_jets._kt2, (JetReconstruction.pt(merged_jet)^2)^_p)
			merged_jet_index = lastindex(_objects)

			ieta_merged_jet, iphi_merged_jet = get_tile(tiling_setup, flat_jets._eta[merged_jet_index],
				flat_jets._phi[merged_jet_index])
			itile_merged_jet = get_tile_linear_index(tiling_setup, ieta_merged_jet, iphi_merged_jet)

			# Set the _sequence for the two merged jets, which is the merged jet index
			push!(_sequences, [_sequences[index_jetA]; _sequences[index_jetB]; merged_jet_index])
			push!(_sequences[index_jetA], merged_jet_index)
			push!(_sequences[index_jetB], merged_jet_index)


			# Delete jetA and jetB from their tiles
			empty!(tainted_slots)
			push!(tainted_slots, index_tile_jetA)
			push!(tainted_slots, index_tile_jetB)
			if itile_merged_jet == index_tile_jetA._itile
				# Put the new jet into jetA's slot
				insert_jet!(tile_jets[itile_merged_jet], index_tile_jetA._ijet, merged_jet_index, flat_jets, _R2)
				index_tile_merged_jet = TiledSoACoord(itile_merged_jet, index_tile_jetA._ijet)
				# Now zap jetB
				push!(tainted_slots, remove_jet!(tile_jets, index_tile_jetB))
			elseif itile_merged_jet == index_tile_jetB._itile
				# Use jetB's slot
				insert_jet!(tile_jets[itile_merged_jet], index_tile_jetB._ijet, merged_jet_index, flat_jets, _R2)
				index_tile_merged_jet = TiledSoACoord(itile_merged_jet, index_tile_jetB._ijet)
				# Now zap jetA
				push!(tainted_slots, remove_jet!(tile_jets, index_tile_jetA))
			else
                # Merged jet is in a different tile
                add_jet!(tile_jets[itile_merged_jet], merged_jet_index, flat_jets, _R2)
                index_tile_merged_jet = TiledSoACoord(itile_merged_jet, tile_jets[itile_merged_jet]._size)
                # Now zap both A and B
                push!(tainted_slots, remove_jet!(tile_jets, index_tile_jetA))
                push!(tainted_slots, remove_jet!(tile_jets, index_tile_jetB))
			end

			# For our new merged jet, scan for nearest neighbours
			# Remember, this is pair-wise, so it will update all jets in its tile and neighbours
			scan_neighbors!(tile_jets, index_tile_merged_jet, _R2)
		end

        # Now take care of tainted neighbours
		empty!(itouched_tiles)
        push!(itouched_tiles, index_tile_jetA._itile)
		union!(itouched_tiles, tile_jets[index_tile_jetA._itile]._nntiles)
        if (jet_merger && index_tile_jetB._itile != index_tile_jetA._itile)
            push!(itouched_tiles, index_tile_jetB._itile)
            union!(itouched_tiles, tile_jets[index_tile_jetB._itile]._nntiles)
        end

		# Scan over the touched tiles, look for jets whose _nn is tainted
		@inbounds for itouched_tile in itouched_tiles
			tile = tile_jets[itouched_tile]
			@inbounds for ijet in 1:tile._size
				if tile._nn[ijet] in tainted_slots
					tile._nn[ijet] = TiledSoACoord(0, 0)
					tile._nndist[ijet] = _R2
					scan_neighbors!(tile_jets, TiledSoACoord(itouched_tile, ijet), _R2)
				end
			end
		end
	end
    # The sequences return value is a list of all jets that merged to this one
	jets, _sequences
end
