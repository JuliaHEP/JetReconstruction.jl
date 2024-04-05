# Structures used internally by the tiled algorithm

"""Structure analogous to BriefJet, but with the extra information
needed for dealing with tiles"""
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
        t = new(-1, 0., 0., 0., 0., -1, 0, 0)
        t.NN = t.previous = t.next = t
        t
    end
    
    TiledJet(id, eta, phi, kt2, NN_dist,
             jet_index, tile_index, dij_posn,
             NN, previous, next) = new(id, eta, phi, kt2, NN_dist,
                                       jet_index, tile_index, dij_posn,
                                       NN, previous, next)
end


const noTiledJet::TiledJet = TiledJet(Nothing)

isvalid(t::TiledJet) = !(t === noTiledJet)

"""
Move a TiledJet in front of a TiledJet list element
The jet to move can be an isolated jet, a jet from another list or a jet from the same list
"""
insert!(nextjet::TiledJet, jettomove::TiledJet) = begin
    if !isnothing(nextjet)
        nextjet.previous  = jettomove
    end

    jettomove.next = nextjet
    jettomove.previous = nextjet.previous
    nextjet = jettomove
end

"""Detach a TiledJet from its list"""
detach!(jet::TiledJet) = begin
    if !isnothing(jet.previous)
        jet.previous.next = jet.next
    end
    if !isnothing(jet.next)
        jet.next.previous = jet.previous
    end
    jet.next = jet.previous = noTiledJet
end

TiledJet(id) = TiledJet(id, 0., 0., 0., 0.,
                        0, 0, 0,
                        noTiledJet, noTiledJet, noTiledJet)

import Base.copy
copy(j::TiledJet) = TiledJet(j.id, j.eta, j.phi, j.kt2, j.NN_dist, j.jets_index, j.tile_index, j.dij_posn, j.NN, j.previous, j.next)

# Iterator over a TiledJet walks along the chain of linked jets
# until we reach a "noTiledJet" (which is !isvalid)
Base.iterate(tj::TiledJet) = begin
   isvalid(tj) ? (tj, tj) : nothing
end
Base.iterate(tj::TiledJet, state::TiledJet) = begin
    isvalid(state.next) ? (state.next::TiledJet, state.next::TiledJet) : nothing
end


"""
Structure with the tiling parameters, as well as some bookkeeping
variables used during reconstruction
"""
struct Tiling
    setup::TilingDef
    tiles::Matrix{TiledJet}
    positions::Matrix{Int}
    tags::Matrix{Bool}
end

"""Return a tiling setup with bookkeeping"""
Tiling(setup::TilingDef) = begin
    t = Tiling(setup,
               fill(noTiledJet, (setup._n_tiles_eta, setup._n_tiles_phi)),
               fill(0, (setup._n_tiles_eta, setup._n_tiles_phi)),
               fill(false, (setup._n_tiles_eta, setup._n_tiles_phi)))
    @inbounds for iphi = 1:setup._n_tiles_phi
        # The order of the following two statements is important
        # to have position = tile_right in case n_tiles_eta = 1
        t.positions[1, iphi] = tile_left
        t.positions[setup._n_tiles_eta, iphi] = tile_right
    end
    t
end

const tile_left = -1
const tile_central = 0
const tile_right = 1

const _n_tile_center = 1
const _n_tile_left_neighbours = 4
const _tile_right_neigbour_indices = 6:9
const _n_tile_right_neighbours = 4
const _n_tile_neighbours = 9

const neigh_init = fill(nothing, _n_tile_neighbours)

"""Structure used for iterating over neighbour tiles"""
struct Surrounding{N}
    indices::NTuple{N, Int}
end

import Base.iterate

Base.iterate(x::T) where {T<:Surrounding} = (x.indices[1], 2)
Base.iterate(x::Surrounding{0}) = nothing
Base.iterate(x::Surrounding{1}, state) = nothing
Base.iterate(x::Surrounding{2}, state) = nothing
Base.iterate(x::Surrounding{3}, state) = state > 3 ? nothing : (x.indices[state], state+1)
Base.iterate(x::Surrounding{4}, state) = state > 4 ? nothing : (x.indices[state], state+1)
Base.iterate(x::Surrounding{6}, state) = state > 6 ? nothing : (x.indices[state], state+1)
Base.iterate(x::Surrounding{9}, state) = state > 9 ? nothing : (x.indices[state], state+1)

import Base.length
length(x::Surrounding{N}) where N = N

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

"""Remove a jet from a tiling"""
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
