# Structures used internally by the tiled algorithm

"""Structure analogous to PseudoJet, but with the extra information
needed for dealing with tiles for the tiled stragegy"""

"""
    struct TiledJet

TiledJet represents a jet in a tiled algorithm for jet reconstruction, with
additional information to track the jet's position in the tiled structures.

# Fields
- `id::Int`: The ID of the jet.
- `eta::Float64`: The rapidity of the jet.
- `phi::Float64`: The azimuthal angle of the jet.
- `kt2::Float64`: The transverse momentum squared of the jet.
- `NN_dist::Float64`: The distance to the nearest neighbor.
- `jets_index::Int`: The index of the jet in the jet array.
- `tile_index::Int`: The index of the tile in the tile array.
- `dij_posn::Int`: The position of this jet in the dij compact array.
- `NN::TiledJet`: The nearest neighbor.
- `previous::TiledJet`: The previous jet.
- `next::TiledJet`: The next jet.
"""
mutable struct TiledJet
    id::Int
    eta::Float64
    phi::Float64
    kt2::Float64
    NN_dist::Float64

    jets_index::Int
    tile_index::Int
    dij_posn::Int

    "Nearest neighbour"
    NN::TiledJet

    previous::TiledJet
    next::TiledJet

    TiledJet(::Type{Nothing}) = begin
        t = new(-1, 0.0, 0.0, 0.0, 0.0, -1, 0, 0)
        t.NN = t.previous = t.next = t
        t
    end

    """
    TiledJet constructor with all fields.
    """
    function TiledJet(id, eta, phi, kt2, NN_dist,
                      jet_index, tile_index, dij_posn,
                      NN, previous, next)
        new(id, eta, phi, kt2, NN_dist,
            jet_index, tile_index, dij_posn,
            NN, previous, next)
    end
end

"""
    const noTiledJet::TiledJet = TiledJet(Nothing)

A constant variable representing a "blank" `TiledJet` object with invalid values.
"""
const noTiledJet::TiledJet = TiledJet(Nothing)

"""
    isvalid(t::TiledJet)

Check if a `TiledJet` is valid, by seeing if it is not the `noTiledJet` object.

# Arguments
- `t::TiledJet`: The `TiledJet` object to check.

# Returns
- `Bool`: `true` if the `TiledJet` object is valid, `false` otherwise.
"""
isvalid(t::TiledJet) = !(t === noTiledJet)

"""
    TiledJet(id)

Constructs a `TiledJet` object with the given `id` and initializes its properties to zero.

# Arguments
- `id`: The ID of the `TiledJet` object.

# Returns
A `TiledJet` object with the specified `id` and values set to zero or noTiledJet.
"""
TiledJet(id) = TiledJet(id, 0.0, 0.0, 0.0, 0.0,
                        0, 0, 0,
                        noTiledJet, noTiledJet, noTiledJet)

"""
    insert!(nextjet::TiledJet, jettomove::TiledJet)

Inserts a `TiledJet` object into the linked list of `TiledJet` objects, before the `nextjet` object.
The jet to move can be an isolated jet, a jet from another list or a jet from the same list

# Arguments
- `nextjet::TiledJet`: The `TiledJet` object after which `jettomove` should be inserted.
- `jettomove::TiledJet`: The `TiledJet` object to be inserted.

# Example
"""
insert!(nextjet::TiledJet, jettomove::TiledJet) = begin
    if !isnothing(nextjet)
        nextjet.previous = jettomove
    end

    jettomove.next = nextjet
    jettomove.previous = nextjet.previous
    nextjet = jettomove
end

"""
    detach!(jet::TiledJet)

Detach a `TiledJet` from its linked list by updating the `previous` and `next` pointers.

# Arguments
- `jet::TiledJet`: The `TiledJet` object to detach.
"""
detach!(jet::TiledJet) = begin
    if !isnothing(jet.previous)
        jet.previous.next = jet.next
    end
    if !isnothing(jet.next)
        jet.next.previous = jet.previous
    end
    jet.next = jet.previous = noTiledJet
end

import Base.copy
"""
    copy(j::TiledJet)

Create a copy of a `TiledJet` object.

# Arguments
- `j::TiledJet`: The `TiledJet` object to be copied.

# Returns
A new `TiledJet` object with the same attributes as the input object.
"""
copy(j::TiledJet) = TiledJet(j.id, j.eta, j.phi, j.kt2, j.NN_dist, j.jets_index,
                             j.tile_index, j.dij_posn, j.NN, j.previous, j.next)

# Iterator over a TiledJet walks along the chain of linked jets
# until we reach a "noTiledJet" (which is !isvalid)

"""
    Base.iterate(tj::TiledJet)

Iterate over a `TiledJet` object's linked list, walking over all jets
until the end (then the next jet is invalid).

# Arguments
- `tj::TiledJet`: The `TiledJet` object to start to iterate over.
"""
Base.iterate(tj::TiledJet) = begin
    isvalid(tj) ? (tj, tj) : nothing
end
function Base.iterate(tj::TiledJet, state::TiledJet)
    isvalid(state.next) ? (state.next::TiledJet, state.next::TiledJet) : nothing
end

"""
    struct Tiling

The `Tiling` struct represents a tiling configuration for jet reconstruction.

# Fields
- `setup::TilingDef`: The tiling definition used for the configuration.
- `tiles::Matrix{TiledJet}`: A matrix of tiled jets, containing the first jet in
  each tile (then the linked list of the first jet is followed to get access to
  all jets in this tile).
- `positions::Matrix{Int}`: Used to track tiles that are on the edge of ϕ array,
  where neighbours need to be wrapped around.
- `tags::Matrix{Bool}`: The matrix of tags indicating whether a tile is valid or
  not (set to `false` initially, then `true` when the tile has been setup
  properly).

"""
struct Tiling
    setup::TilingDef
    tiles::Matrix{TiledJet}
    positions::Matrix{Int}
    tags::Matrix{Bool}
end

"""
    Tiling(setup::TilingDef)

Constructs a intial `Tiling` object based on the provided `setup` parameters.

# Arguments
- `setup::TilingDef`: The setup parameters for the tiling.

# Returns
A `Tiling` object.
"""
Tiling(setup::TilingDef) = begin
    t = Tiling(setup,
               fill(noTiledJet, (setup._n_tiles_eta, setup._n_tiles_phi)),
               fill(0, (setup._n_tiles_eta, setup._n_tiles_phi)),
               fill(false, (setup._n_tiles_eta, setup._n_tiles_phi)))
    @inbounds for iphi in 1:(setup._n_tiles_phi)
        # The order of the following two statements is important
        # to have position = tile_right in case n_tiles_eta = 1
        t.positions[1, iphi] = tile_left
        t.positions[setup._n_tiles_eta, iphi] = tile_right
    end
    t
end

"Signal that a tile is on the left hand side of the tiling array in ϕ"
const tile_left = -1
"Signal that a tile is central for the tiling array in ϕ"
const tile_central = 0
"Signal that a tile is on the right hand side of the tiling array in ϕ"
const tile_right = 1
"Number of neighbours for a tile in the tiling array (including itself)"
const _n_tile_neighbours = 9

"""
    struct Surrounding{N}

Structure used for iterating over neighbour tiles.

# Fields
- `indices::NTuple{N, Int}`: A tuple of `N` integers representing the indices.
"""
struct Surrounding{N}
    indices::NTuple{N, Int}
end

import Base.iterate

Base.iterate(x::T) where {T <: Surrounding} = (x.indices[1], 2)
Base.iterate(x::Surrounding{0}) = nothing
Base.iterate(x::Surrounding{1}, state) = nothing
Base.iterate(x::Surrounding{2}, state) = nothing
Base.iterate(x::Surrounding{3}, state) = state > 3 ? nothing : (x.indices[state], state + 1)
Base.iterate(x::Surrounding{4}, state) = state > 4 ? nothing : (x.indices[state], state + 1)
Base.iterate(x::Surrounding{6}, state) = state > 6 ? nothing : (x.indices[state], state + 1)
Base.iterate(x::Surrounding{9}, state) = state > 9 ? nothing : (x.indices[state], state + 1)

import Base.length
length(x::Surrounding{N}) where {N} = N

"""
    surrounding(center::Int, tiling::Tiling)

Compute the surrounding indices of a given center index in a tiling.

# Arguments
- `center::Int`: The center index.
- `tiling::Tiling`: The tiling object.

# Returns
- `Surrounding`: An object containing the surrounding indices.
"""
surrounding(center::Int, tiling::Tiling) = begin
    #                        4|6|9
    #                        3|1|8
    #                        2|5|7
    #  -> η

    iphip = mod1(center + tiling.setup._n_tiles_eta, tiling.setup._n_tiles)
    iphim = mod1(center - tiling.setup._n_tiles_eta, tiling.setup._n_tiles)

    if tiling.setup._n_tiles_eta == 1
        return Surrounding{3}((center, iphim, iphip))
    elseif tiling.positions[center] == tile_right
        return Surrounding{6}((center, iphim, iphip, iphim - 1, center - 1, iphip - 1))
    elseif tiling.positions[center] == tile_central
        return Surrounding{9}((center, iphim - 1, center - 1, iphip - 1,
                               iphim, iphip,
                               iphim + 1, center + 1, iphip + 1))
    else #tile_left
        return Surrounding{6}((center, iphim, iphip,
                               iphim + 1, center + 1, iphip + 1))
    end
end

"""
    rightneighbours(center::Int, tiling::Tiling)

Compute the indices of the right neighbors of a given center index in a tiling.
This is used in the inital sweep to calculate the nearest neighbors, where the
search between jets for the nearest neighbour is bi-directional, thus when a
tile is considered only the right neighbours are needed to compare jet
distances as the left-hand tiles have been done from that tile already.

# Arguments
- `center::Int`: The center index.
- `tiling::Tiling`: The tiling object.

# Returns
- `Surrounding`: An object containing the indices of the right neighbors.
"""
rightneighbours(center::Int, tiling::Tiling) = begin
    #                         |1|4
    #                         | |3
    #                         | |2
    #  -> η

    iphip = mod1(center + tiling.setup._n_tiles_eta, tiling.setup._n_tiles)
    iphim = mod1(center - tiling.setup._n_tiles_eta, tiling.setup._n_tiles)

    if tiling.positions[center] == tile_right
        return Surrounding{1}((iphip,))
    else
        return Surrounding{4}((iphip, iphim + 1, center + 1, iphip + 1))
    end
end

"""
    tiledjet_remove_from_tiles!(tiling, jet)

Remove a jet from the given tiling structure.

# Arguments
- `tiling`: The tiling structure from which the jet will be removed.
- `jet`: The jet to be removed from the tiling structure.

# Description
This function removes a jet from the tiling structure. It adjusts the linked list
to be consistent with the removal of the jet.
"""
tiledjet_remove_from_tiles!(tiling, jet) = begin
    if !isvalid(jet.previous)
        # We are at head of the tile, so reset it.
        # If this was the only jet on the tile then tile->head will now be NULL
        tiling.tiles[jet.tile_index] = jet.next
    else
        # Adjust link from previous jet in this tile
        jet.previous.next = jet.next
    end

    if isvalid(jet.next)
        # Adjust backwards-link from next jet in this tile
        jet.next.previous = jet.previous
    end

    jet.next = jet.previous = noTiledJet # To be clean, but not mandatory
end
