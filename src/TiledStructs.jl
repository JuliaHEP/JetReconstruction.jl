# Defined the structures and associated functions used in tiled
# jet reconstruction

import Base.==
import Base.hash

"""Tiling definition parameters"""
struct TilingDef
	_tiles_eta_min::Float64  # Minimum rapidity
	_tiles_eta_max::Float64  # Maximum rapidity
	_tile_size_eta::Float64  # Size of a tile in rapidity (usually R^2)
	_tile_size_phi::Float64  # Size of a tile in phi (usually a bit more than R^2)
	_n_tiles_eta::Int   # Number of tiles across rapidity
	_n_tiles_phi::Int   # Number of tiles across phi
	_n_tiles::Int       # Total number of tiles
	_tiles_ieta_min::Int # Min_rapidity / rapidity tile size (needed?)
	_tiles_ieta_max::Int # Max_rapidity / rapidity tile size (needed?)

	# This maps (ieta, iphi) to the linear index of the tile jets
	# _tile_linear_indexes

	# And back again...
	# _tile_cartesian_indexes
	
	# Use an inner constructor as _n_tiles and _tile_linear_indexes 
	# are defined by the other values
	function TilingDef(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
		_n_tiles_eta, _n_tiles_phi, _tiles_ieta_min, _tiles_ieta_max)
		new(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
		_n_tiles_eta, _n_tiles_phi, _n_tiles_eta*_n_tiles_phi, _tiles_ieta_min, _tiles_ieta_max,
		# LinearIndices((1:_n_tiles_eta, 1:_n_tiles_phi)), CartesianIndices((1:_n_tiles_eta, 1:_n_tiles_phi))
		)
	end
end

"""Nearest neighbour coordinates"""
mutable struct TiledSoACoord
    _itile::Int           # Jet tile index (flattened)
    _ijet::Int            # Jet position in this tile
end

"""Correct hash for the coordinate pair"""
hash(x::TiledSoACoord, h::UInt) = begin
	my_hash = hash(hash(x._itile, UInt(0)), hash(x._ijet, UInt(0)))
	hash(my_hash, h)
end

"""Equality operator for tiled coordinates"""
==(x::TiledSoACoord,y::TiledSoACoord) = hash(x, UInt(0))==hash(y, UInt(0))

"""Setter for nearest neighbour"""
set_nn!(mynn::TiledSoACoord, itile, ijet) = begin
    mynn._itile = itile
    mynn._ijet = ijet
end

"""Do we have a NN, or not"""
valid_nn(mynn::TiledSoACoord) = begin
    return mynn._itile > 0
end

abstract type JetSoA end

"""Structure of arrays for tiled jet parameters, using an SoA layout
for computational efficiency"""
mutable struct TiledJetSoA <: JetSoA
    _size::Int              # Active jet count (can be less than the vector length)
	_kt2::Vector{Float64}         # p_t^-2p
	_eta::Vector{Float64}         # Rapidity
	_phi::Vector{Float64}         # Phi coordinate
	_index::Vector{Int}     # My jet index
	_nn::Vector{TiledSoACoord}    # Nearest neighbour location (if (0,0) no nearest neighbour)
	_nndist::Vector{Float64}      # Distance to my nearest neighbour
    _dij::Vector{Float64}         # Jet metric distance to my nearest neighbour
    _righttiles::Vector{Int} # Indexes of all tiles to my right
    _nntiles::Vector{Int}   # Indexes of all neighbour tiles
	# Inner constructor to ensure the sizehint for tile caches
	function TiledJetSoA(n::Int)
		my_tiledjet = new(n,
			Vector{Float64}(undef, n),
			Vector{Float64}(undef, n),
			Vector{Float64}(undef, n),
			Vector{Int}(undef, n),
			Vector{TiledSoACoord}(undef, n),
			Vector{Float64}(undef, n),
			Vector{Float64}(undef, n),
			Vector{Int}(undef, 0),
			Vector{Int}(undef, 0)
		)
		sizehint!(my_tiledjet._righttiles, 4)
		sizehint!(my_tiledjet._nntiles, 8)
		my_tiledjet
	end
end

"""Return the flat jet index of the nearest neighbour tile of a jet"""
nnindex(tile_jets::Array{TiledJetSoA, 2}, itile, ijet) = begin
	if tile_jets[itile]._nn[ijet]._itile == 0
		return 0
	end
    return tile_jets[tile_jets[itile]._nn[ijet]._itile]._index[tile_jets[itile]._nn[ijet]._ijet]
end

"""Structure for the flat jet SoA, as it's convenient"""
mutable struct FlatJetSoA <: JetSoA
	_size::Int            # Number of active entries (may be less than the vector size!)
	_kt2::Vector{Float64}       # p_t^-2
	_eta::Vector{Float64}       # Rapidity
	_phi::Vector{Float64}       # Phi coordinate
	_index::Vector{Int}   # My jet index
end

"""Insert a jet into a tile at the given slot"""
insert_jet!(tile::TiledJetSoA, slot::Int, index::Int, flat_jets::FlatJetSoA, R2::AbstractFloat) = begin
	tile._kt2[slot] = flat_jets._kt2[index]
	tile._eta[slot] = flat_jets._eta[index]
	tile._phi[slot] = flat_jets._phi[index]
	tile._index[slot] = index
	tile._nn[slot] = TiledSoACoord(0, 0)
	tile._nndist[slot] = R2
	tile._dij[slot] = R2 * tile._kt2[slot]
end

"""Add a jet to a tile beyond the active slots"""
add_jet!(tile::TiledJetSoA, index::Int, flat_jets::FlatJetSoA, R2::AbstractFloat) = begin
	tile._size += 1
	push!(tile._kt2, flat_jets._kt2[index])
	push!(tile._eta, flat_jets._eta[index])
	push!(tile._phi, flat_jets._phi[index])
	push!(tile._index, index)
	push!(tile._nn, TiledSoACoord(0, 0))
	push!(tile._nndist, R2)
	push!(tile._dij, R2 * flat_jets._kt2[index])
end

"""Remove a jet from the tile at the given tile index and repack if needed"""
remove_jet!(tile_jets::Array{TiledJetSoA, 2}, tile_index::TiledSoACoord) = begin
	tile = tile_jets[tile_index._itile]
	if tile._size != tile_index._ijet
		# Need to copy the last jet into the slot, to ensure a contiguous array
		# These two slots become tainted
		tile._kt2[tile_index._ijet] = tile._kt2[tile._size]
		tile._eta[tile_index._ijet] = tile._eta[tile._size]
		tile._phi[tile_index._ijet] = tile._phi[tile._size]
		tile._index[tile_index._ijet] = tile._index[tile._size]
		tile._nn[tile_index._ijet] = tile._nn[tile._size]
		tile._nndist[tile_index._ijet] = tile._nndist[tile._size]
		tile._dij[tile_index._ijet] = tile._dij[tile._size]
		tainted = TiledSoACoord(tile_index._itile, tile._size)
	else
		tainted = TiledSoACoord(0, 0)
	end
	# Physically reduce the arrays - safer while developing to avoid accidents!
	# For optimised code, consider not doing this, as we can use _size as a bound when scanning
	deleteat!(tile._kt2, tile._size)
	deleteat!(tile._eta, tile._size)

	deleteat!(tile._phi, tile._size)
	deleteat!(tile._index, tile._size)
	deleteat!(tile._nn, tile._size)
	deleteat!(tile._nndist, tile._size)
	deleteat!(tile._dij, tile._size)
	# Removing the last entry
	tile._size -= 1
	return tainted
end

get_jet(tile_jets::Array{TiledJetSoA, 2}, tile_index::TiledSoACoord) = begin
	tile = tile_jets[tile_index._itile]
	ijet = tile_index._ijet
	return "($(tile._eta[ijet]), $(tile._phi[ijet]) [$(tile._index[ijet])]) -> $(tile._nn[ijet])"
end


"""Return the nth jet in the SoA"""
get_jet(j, n::Int) = begin
    return j._index[n], j._eta[n], j._phi[n], j._kt2[n]
end

# Getters - will work for both SoAs
kt2(j::JetSoA, n::Int) = j._kt2[n]
eta(j::JetSoA, n::Int) = j._eta[n]
phi(j::JetSoA, n::Int) = j._phi[n]
index(j::JetSoA, n::Int) = j._index[n]
nn(j::JetSoA, n::Int) = j._nn[n]
nndist(j::JetSoA, n::Int) = j._nndist[n]

# Setters
set_kt2!(j::JetSoA, n::Int, v) = j._kt2[n] = v
set_eta!(j::JetSoA, n::Int, v) = j._eta[n] = v
set_phi!(j::JetSoA, n::Int, v) = j._phi[n] = v
set_index!(j::JetSoA, n::Int, v) = j._index[n] = v
set_nn!(j::JetSoA, n::Int, v::TiledSoACoord) = j._nn[n] = v
set_nndist!(j::JetSoA, n::Int, v) = j._nndist[n] = v
