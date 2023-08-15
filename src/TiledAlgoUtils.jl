# Common functions and structures for Tiled reconstruction algorithms

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
	
	# Use an inner constructor as _n_tiles and _tile_linear_indexes 
	# are defined by the other values
	function TilingDef(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
		_n_tiles_eta, _n_tiles_phi, _tiles_ieta_min, _tiles_ieta_max)
		new(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
		_n_tiles_eta, _n_tiles_phi, _n_tiles_eta*_n_tiles_phi, _tiles_ieta_min, _tiles_ieta_max,
		)
	end
end

"""
Determine an extent for the rapidity tiles, by binning in 
rapidity and collapsing the outer bins until they have about
1/4 the number of particles as the maximum bin. This is the 
heuristic which is used by FastJet.
"""
function determine_rapidity_extent(_eta::Vector{T}) where T <: AbstractFloat
	length(_eta) == 0 && return 0.0, 0.0

	nrap = 20
	nbins = 2 * nrap
	counts = zeros(Int, nbins)

	# Get the minimum and maximum rapidities and at the same time bin
	# the multiplicities as a function of rapidity to help decide how
	# far out it's worth going
	minrap = floatmax(T)
	maxrap = -floatmax(T)
	ibin = 0
	for y in _eta
		minrap = min(minrap, y)
		maxrap = max(maxrap, y)

		# Bin in rapidity to decide how far to go with the tiling.
		# The bins go from ibin=1 (rap=-infinity..-19)
		# to ibin = nbins (rap=19..infinity for nrap=20)
		ibin = clamp(1 + unsafe_trunc(Int, y + nrap), 1, nbins)
		@inbounds counts[ibin] += 1
	end

	# Get the busiest bin
	max_in_bin = maximum(counts)

	# Now find minrap, maxrap such that edge bin never contains more
	# than some fraction of busiest, and at least a few particles; first do
	# it from left. NB: the thresholds chosen here are largely
	# guesstimates as to what might work.
	#
	# 2014-07-17: in some tests at high multiplicity (100k) and particles going up to
	#             about 7.3, anti-kt R=0.4, we found that 0.25 gave 20% better run times
	#             than the original value of 0.5.
	allowed_max_fraction = 0.25

	# The edge bins should also contain at least min_multiplicity particles
	min_multiplicity = 4

	# now calculate how much we can accumulate into an edge bin
	allowed_max_cumul = floor(max(max_in_bin * allowed_max_fraction,
		min_multiplicity))

	# make sure we don't require more particles in a bin than max_in_bin
	allowed_max_cumul = min(max_in_bin, allowed_max_cumul)

	# start scan over rapidity bins from the "left", to find out minimum rapidity of tiling
	cumul_lo = 0.0
	ibin_lo = 1
	while ibin_lo <= nbins
		@inbounds cumul_lo += counts[ibin_lo]
		if cumul_lo >= allowed_max_cumul
			minrap = max(minrap, ibin_lo - nrap - 1)
			break
		end
		ibin_lo += 1
	end
	@assert ibin_lo != nbins # internal consistency check that you found a bin

	# then do it from "right", to find out maximum rapidity of tiling
	cumul_hi = 0.0
	ibin_hi = nbins
	while ibin_hi >= 1
		@inbounds cumul_hi += counts[ibin_hi]
		if cumul_hi >= allowed_max_cumul
			maxrap = min(maxrap, ibin_hi - nrap)
			break
		end
		ibin_hi -= 1
	end
	@assert ibin_hi >= 1 # internal consistency check that you found a bin

	# consistency check
	@assert ibin_hi >= ibin_lo

	minrap, maxrap
end

"""
Setup the tiling parameters for this recontstruction
"""
function setup_tiling(eta::Vector{T}, Rparam::AbstractFloat) where T <: AbstractFloat
	# First decide tile sizes (with a lower bound to avoid huge memory use with
	# very small R)
	tile_size_eta = max(0.1, Rparam)

	# It makes no sense to go below 3 tiles in phi -- 3 tiles is
	# sufficient to make sure all pair-wise combinations up to pi in
	# phi are possible
	n_tiles_phi   = max(3, floor(Int, 2π / tile_size_eta))
	tile_size_phi = 2π / n_tiles_phi # >= Rparam and fits in 2pi

	tiles_eta_min, tiles_eta_max = determine_rapidity_extent(eta)

	# now adjust the values
	tiles_ieta_min = floor(Int, tiles_eta_min / tile_size_eta)
	tiles_ieta_max = floor(Int, tiles_eta_max / tile_size_eta) #FIXME shouldn't it be ceil ?
	tiles_eta_min = tiles_ieta_min * tile_size_eta
	tiles_eta_max = tiles_ieta_max * tile_size_eta
	n_tiles_eta = tiles_ieta_max - tiles_ieta_min + 1

	tiling_setup = TilingDef(tiles_eta_min, tiles_eta_max,
		tile_size_eta, tile_size_phi,
		n_tiles_eta, n_tiles_phi,
		tiles_ieta_min, tiles_ieta_max)

	tiling_setup
end

"""
Return the geometric distance between a pair of (eta,phi) coordinates
"""
geometric_distance(eta1::AbstractFloat, phi1::AbstractFloat, eta2::AbstractFloat, phi2::AbstractFloat) = begin
	δeta = eta2 - eta1
	δphi = π - abs(π - abs(phi1 - phi2))
	return δeta * δeta + δphi * δphi
end

"""
Return the dij metric distance between a pair of pseudojets
"""
get_dij_dist(nn_dist, kt2_1, kt2_2, R2) = begin
	if kt2_2 == 0.0
		return kt2_1 * R2
	end
	return nn_dist * min(kt2_1, kt2_2)
end

"""
Return the tile coordinates of an (eta, phi) pair
"""
get_tile(tiling_setup::TilingDef, eta::AbstractFloat, phi::AbstractFloat) = begin
	# The eta clamp is necessary as the extreme bins catch overflows for high abs(eta)
	ieta = clamp(floor(Int, (eta - tiling_setup._tiles_eta_min) / tiling_setup._tile_size_eta), 1, tiling_setup._n_tiles_eta)
	# The phi clamp should not really be necessary, as long as phi values are [0,2π)
	iphi = clamp(floor(Int, 1 + (phi / 2π) * tiling_setup._n_tiles_phi), 1, tiling_setup._n_tiles_phi)
	ieta, iphi
end

"""
Map an (η, ϕ) pair into a linear index, which is much faster "by hand" than using
the LinearIndices construct (like x100, which is bonkers, but there you go...)
"""
get_tile_linear_index(tiling_setup::TilingDef, i_η::Int, i_ϕ::Int) = begin
	return tiling_setup._n_tiles_eta * (i_ϕ-1) + i_η
end

"""
Map a linear index to a tuple of (η, ϕ) - again this is very much faster than
using CartesianIndices
"""
get_tile_cartesian_indices(tiling_setup::TilingDef, index::Int) = begin
	return (rem(index-1, tiling_setup._n_tiles_eta)+1, div(index-1, tiling_setup._n_tiles_eta)+1)
end

"""
Iterator for the indexes of rightmost tiles for a given Cartesian tile index
	- These are the tiles above and to the right of the given tile (X=yes, O=no)
		XXX
		O.X
		OOO
	- η coordinate must be in range, ϕ coordinate wraps

"""
struct rightmost_tiles
    n_η::Int		# Number of η tiles
    n_ϕ::Int		# Number of ϕ tiles
    start_η::Int	# Centre η tile coordinate
    start_ϕ::Int	# Centre ϕ tile coordinate
end

function Base.iterate(t::rightmost_tiles, state=1)
    mapping = ((-1,-1), (-1,0), (-1,1), (0,1))
    if t.start_η == 1 && state == 1
        state = 4
    end
    while state <= 4
        η = t.start_η + mapping[state][1]
        ϕ = t.start_ϕ + mapping[state][2]
        if ϕ > t.n_ϕ
            ϕ = 1
        elseif ϕ < 1
            ϕ = t.n_ϕ
        end
        return (η,ϕ), state+1
    end
    return nothing
end

"""
Iterator for the indexes of neighbouring tiles for a given Cartesian tile index
		XXX
		X.X
		XXX
	- η coordinate must be in range, ϕ coordinate wraps

"""
struct neighbour_tiles
    n_η::Int		# Number of η tiles
    n_ϕ::Int		# Number of ϕ tiles
    start_η::Int	# Centre η tile coordinate
    start_ϕ::Int	# Centre ϕ tile coordinate
end

function Base.iterate(t::neighbour_tiles, state=1)
    mapping = ((-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1))
	# Skip for top row in η
    if t.start_η == 1 && state == 1
        state = 4
    end
	# Skip for bottom row in η
	if t.start_η == t.n_η && state == 6
		state = 9
	end
    while state <= 8
        η = t.start_η + mapping[state][1]
        ϕ = t.start_ϕ + mapping[state][2]
        if ϕ > t.n_ϕ
            ϕ = 1
        elseif ϕ < 1
            ϕ = t.n_ϕ
        end
        return (η,ϕ), state+1
    end
    return nothing
end
