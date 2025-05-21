# Common functions and structures for Tiled reconstruction algorithms

"""
	struct TilingDef

A struct representing the definition of a specific tiling scheme.

# Fields
- `_tiles_eta_min::Float64`: The minimum rapidity of the tiles.
- `_tiles_eta_max::Float64`: The maximum rapidity of the tiles.
- `_tile_size_eta::Float64`: The size of a tile in rapidity (usually R^2).
- `_tile_size_phi::Float64`: The size of a tile in phi (usually a bit more than R^2).
- `_n_tiles_eta::Int`: The number of tiles across rapidity.
- `_n_tiles_phi::Int`: The number of tiles across phi.
- `_n_tiles::Int`: The total number of tiles.
- `_tiles_ieta_min::Int`: The minimum rapidity tile index.
- `_tiles_ieta_max::Int`: The maximum rapidity tile index.

# Constructor
	TilingDef(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
		_n_tiles_eta, _n_tiles_phi, _tiles_ieta_min, _tiles_ieta_max)

Constructs a `TilingDef` object with the given parameters.
"""
struct TilingDef
    _tiles_eta_min::Float64
    _tiles_eta_max::Float64
    _tile_size_eta::Float64
    _tile_size_phi::Float64
    _n_tiles_eta::Int
    _n_tiles_phi::Int
    _n_tiles::Int
    _tiles_ieta_min::Int
    _tiles_ieta_max::Int

    function TilingDef(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
                       _n_tiles_eta, _n_tiles_phi, _tiles_ieta_min, _tiles_ieta_max)
        new(_tiles_eta_min, _tiles_eta_max, _tile_size_eta, _tile_size_phi,
            _n_tiles_eta, _n_tiles_phi, _n_tiles_eta * _n_tiles_phi, _tiles_ieta_min,
            _tiles_ieta_max)
    end
end

"""
    determine_rapidity_extent(eta::Vector{T}) where T <: AbstractFloat

Calculate the minimum and maximum rapidities based on the input vector `eta`.
The function determines the rapidity extent by binning the multiplicities as a
function of rapidity and finding the minimum and maximum rapidities such that
the edge bins contain a certain fraction (~1/4) of the busiest bin and a minimum
number of particles.

This is the heuristic which is used by FastJet (inline comments are from FastJet).

## Arguments
- `eta::Vector{T}`: A vector of rapidity values.

## Returns
- `minrap::T`: The minimum rapidity value.
- `maxrap::T`: The maximum rapidity value.
"""
function determine_rapidity_extent(eta::Vector{T}) where {T <: AbstractFloat}
    length(eta) == 0 && return 0.0, 0.0

    nrap = 20
    nbins = 2 * nrap
    counts = zeros(Int, nbins)

    # Get the minimum and maximum rapidities and at the same time bin
    # the multiplicities as a function of rapidity to help decide how
    # far out it is worth going
    minrap = maxrap = eta[1]
    ibin = 0
    for y in eta
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
    setup_tiling(eta::Vector{T}, Rparam::AbstractFloat) where T <: AbstractFloat

This function sets up the tiling parameters for a reconstruction given a vector
of rapidities `eta` and a radius parameter `Rparam`.

# Arguments
- `eta::Vector{T}`: A vector of rapidities.
- `Rparam::AbstractFloat`: The jet radius parameter.

# Returns
- `tiling_setup`: A `TilingDef` object containing the tiling setup parameters.

# Description
The function first decides the tile sizes based on the `Rparam` value. It then
determines the number of tiles in the phi direction (`n_tiles_phi`) based on the
tile size. Next, it determines the rapidity extent of the input `eta` vector and
adjusts the values accordingly. Finally, it creates a `TilingDef` object with
the calculated tiling parameters and returns it.
"""
function setup_tiling(eta::Vector{T}, Rparam::AbstractFloat) where {T <: AbstractFloat}
    # First decide tile sizes (with a lower bound to avoid huge memory use with
    # very small R)
    tile_size_eta = max(0.1, Rparam)

    # It makes no sense to go below 3 tiles in phi -- 3 tiles is
    # sufficient to make sure all pair-wise combinations up to pi in
    # phi are possible
    n_tiles_phi = max(3, floor(Int, 2π / tile_size_eta))
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
	geometric_distance(eta1::AbstractFloat, phi1::AbstractFloat, eta2::AbstractFloat, phi2::AbstractFloat)

Compute the geometric distance between two points in the rap-phi plane.

# Arguments
- `eta1::AbstractFloat`: The eta coordinate of the first point.
- `phi1::AbstractFloat`: The phi coordinate of the first point.
- `eta2::AbstractFloat`: The eta coordinate of the second point.
- `phi2::AbstractFloat`: The phi coordinate of the second point.

# Returns
- `distance::Float64`: The geometric distance between the two points.
"""
function geometric_distance(eta1::AbstractFloat, phi1::AbstractFloat, eta2::AbstractFloat,
                            phi2::AbstractFloat)
    δeta = eta2 - eta1
    δphi = π - abs(π - abs(phi1 - phi2))
    return δeta * δeta + δphi * δphi
end

"""
	get_dij_dist(nn_dist, kt2_1, kt2_2, R2)

Compute the dij metric distance between two jets.

# Arguments
- `nn_dist`: The nearest-neighbor distance between two jets.
- `kt2_1`: The squared momentum metric value of the first jet.
- `kt2_2`: The squared momentum metric value of the second jet.
- `R2`: The jet radius parameter squared.

# Returns
The distance between the two jets.

If `kt2_2` is equal to 0.0, then the first jet doesn't actually have a valid 
neighbour, so it's treated as a single jet adjacent to the beam.
"""
get_dij_dist(nn_dist, kt2_1, kt2_2, R2) = begin
    if kt2_2 == 0.0
        return kt2_1 * R2
    end
    return nn_dist * min(kt2_1, kt2_2)
end

"""
    get_tile(tiling_setup::TilingDef, eta::AbstractFloat, phi::AbstractFloat)

Given a `tiling_setup` object, `eta` and `phi` values, this function calculates
the tile indices for the given `eta` and `phi` values.

# Arguments
- `tiling_setup`: A `TilingDef` object that contains the tiling setup
  parameters.
- `eta`: The eta value for which to calculate the tile index.
- `phi`: The phi value for which to calculate the tile index.

# Returns
- `ieta`: The tile index along the eta direction.
- `iphi`: The tile index along the phi direction.
"""
function get_tile(tiling_setup::TilingDef, eta::AbstractFloat, phi::AbstractFloat)
    # The eta clamp is necessary as the extreme bins catch overflows for high abs(eta)
    ieta = clamp(floor(Int,
                       (eta - tiling_setup._tiles_eta_min) / tiling_setup._tile_size_eta),
                 1, tiling_setup._n_tiles_eta)
    # The phi clamp should not really be necessary, as long as phi values are [0,2π)
    iphi = clamp(floor(Int, 1 + (phi / 2π) * tiling_setup._n_tiles_phi), 1,
                 tiling_setup._n_tiles_phi)
    ieta, iphi
end

"""
    get_tile_linear_index(tiling_setup::TilingDef, i_η::Int, i_ϕ::Int)

Compute the linear index of a tile in a tiled setup. This is much faster in this
function than using the LinearIndices construct (like x100, which is bonkers,
but there you go...)

# Arguments
- `tiling_setup::TilingDef`: The tiling setup defining the number of tiles in
  each dimension.
- `i_η::Int`: The index of the tile in the η dimension.
- `i_ϕ::Int`: The index of the tile in the ϕ dimension.

# Returns
- The linear index of the tile.
"""
function get_tile_cartesian_indices(tiling_setup::TilingDef, index::Int)
    return (rem(index - 1, tiling_setup._n_tiles_eta) + 1,
            div(index - 1, tiling_setup._n_tiles_eta) + 1)
end

"""
Iterator for the indexes of rightmost tiles for a given Cartesian tile index
	- These are the tiles above and to the right of the given tile (X=yes, O=no)
		XXX
		O.X
		OOO
	- rapidity coordinate must be in range, ϕ coordinate wraps

"""

"""
    struct rightmost_tiles

A struct for iterating over rightmost tiles for a given Cartesian tile index.
These are the tiles above and to the right of the given tile (X=included, O=not
included):

	XXX
	O.X
	OOO

Note, rapidity coordinate must be in range, ϕ coordinate wraps

# Fields
- `n_η::Int`: Number of η tiles
- `n_ϕ::Int`: Number of ϕ tiles
- `start_η::Int`: Centre η tile coordinate
- `start_ϕ::Int`: Centre ϕ tile coordinate
"""
struct rightmost_tiles
    n_η::Int# Number of η tiles
    n_ϕ::Int# Number of ϕ tiles
    start_η::Int# Centre η tile coordinate
    start_ϕ::Int# Centre ϕ tile coordinate
end

"""
    Base.iterate(t::rightmost_tiles, state=1)

Iterate over the `rightmost_tiles` object, returning all the rightmost tiles for
a given Cartesian tile index.
"""
function Base.iterate(t::rightmost_tiles, state = 1)
    mapping = ((-1, -1), (-1, 0), (-1, 1), (0, 1))
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
        return (η, ϕ), state + 1
    end
    return nothing
end

"""
    struct neighbour_tiles

A struct representing the neighbouring tiles.

A struct for iterating over all neighbour tiles for a given Cartesian tile index.
These are the tiles above and to the right of the given tile (X=included, O=not
included):

	XXX
	X.X
	XXX

Note, rapidity coordinate must be in range, ϕ coordinate wraps

# Fields
- `n_η::Int`: Number of η tiles
- `n_ϕ::Int`: Number of ϕ tiles
- `start_η::Int`: Centre η tile coordinate
- `start_ϕ::Int`: Centre ϕ tile coordinate
"""
struct neighbour_tiles
    n_η::Int# Number of η tiles
    n_ϕ::Int# Number of ϕ tiles
    start_η::Int# Centre η tile coordinate
    start_ϕ::Int# Centre ϕ tile coordinate
end

"""
    Base.iterate(t::neighbour_tiles, state=1)

Iterate over the `neighbour_tiles` object, returning all the neighbour tiles for
a given Cartesian tile index.
"""
function Base.iterate(t::neighbour_tiles, state = 1)
    mapping = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
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
        return (η, ϕ), state + 1
    end
    return nothing
end
