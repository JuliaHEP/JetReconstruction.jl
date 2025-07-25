using Statistics

"""
    struct SoftKiller

Implements the SoftKiller pileup mitigation algorithm, inspired by the FastJet contrib package.

The original algorithm is described in:
https://arxiv.org/abs/1407.0408
by Matteo Cacciari, Gavin P. Salam, Gregory Soyez

This version inspired by the SoftKiller implementation in the FastJet contrib package
Original C++ code: https://fastjet.hepforge.org/contrib/

The SoftKiller algorithm divides the rapidity-phi plane into a grid of tiles and determines a dynamic
pt threshold by finding the median of the maximum pt in each tile. Particles with pt below this threshold
are removed from the event.

# Fields
- `_ymax::Float64`: Maximum rapidity of the grid.
- `_ymin::Float64`: Minimum rapidity of the grid.
- `_requested_drap::Float64`: Requested grid spacing in rapidity.
- `_requested_dphi::Float64`: Requested grid spacing in phi.
- `_ntotal::Int`: Total number of tiles.
- `_dy::Float64`: Actual grid spacing in rapidity.
- `_dphi::Float64`: Actual grid spacing in phi.
- `_cell_area::Float64`: Area of a single tile.
- `_inverse_dy::Float64`: Inverse of rapidity grid spacing.
- `_inverse_dphi::Float64`: Inverse of phi grid spacing.
- `_ny::Int`: Number of tiles in rapidity.
- `_nphi::Int`: Number of tiles in phi.

# Constructors
- `SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)`: 
  Construct a grid from `rapmin` to `rapmax` in rapidity, with tile sizes `drap` and `dphi`.
- `SoftKiller(rapmax::Float64, grid_size::Float64)`: 
  Construct a square grid from `-rapmax` to `rapmax` in rapidity, with tile size `grid_size`.

"""
struct SoftKiller
    _ymin::Float64
    _ymax::Float64
    _requested_drap::Float64
    _requested_dphi::Float64
    _ntotal::Int
    _dy::Float64
    _dphi::Float64
    _cell_area::Float64
    _inverse_dy::Float64
    _inverse_dphi::Float64
    _ny::Int
    _nphi::Int

    """
        SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)

    Construct a SoftKiller grid from `rapmin` to `rapmax` in rapidity, with tile sizes `drap` and `dphi`.
    """

    function SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)
        @assert rapmax > rapmin
        @assert drap > 0
        @assert dphi > 0

        ny_double = (rapmax - rapmin) / drap
        ny = max(round(Int, ny_double + 0.5), 1)
        dy = (rapmax - rapmin) / ny
        inverse_dy = ny / (rapmax - rapmin)

        nphi = round(Int, (2 * π) / dphi + 0.5)
        dphi_final = (2 * π) / nphi
        inverse_dphi = nphi / (2 * π)

        @assert ny >= 1 && nphi >= 1

        ntotal = nphi * ny
        cell_area = dy * dphi_final

        new(rapmin, rapmax, drap, dphi, ntotal, dy, dphi_final, cell_area,
            inverse_dy, inverse_dphi, ny, nphi)
    end

    """
        SoftKiller(rapmax::Float64, grid_size::Float64)

    Construct a square SoftKiller grid from `-rapmax` to `rapmax` in rapidity, with tile size `grid_size`.
    """
    function SoftKiller(rapmax::Float64, grid_size::Float64)
        return SoftKiller(-rapmax, rapmax, grid_size, grid_size)
    end
end

"""
    tile_index(sk::SoftKiller, p::PseudoJet)

Return the tile index for a given `PseudoJet` `p` in the SoftKiller grid `sk`.
Returns -1 if the jet is outside the grid bounds.
"""
function tile_index(sk::SoftKiller, p::PseudoJet)
    y_minus_ymin = rapidity(p) - sk._ymin
    if y_minus_ymin < 0
        return -1
    end

    iy = round(Int, y_minus_ymin * sk._inverse_dy)
    if iy >= sk._ny
        return -1
    end

    iphi = round(Int, phi(p) * sk._inverse_dphi)
    if iphi == sk._nphi
        iphi = 0
    end

    res = round(Int, iy * sk._nphi + iphi)

    res + 1
end

import Base: show

"""
    show(io::IO, sk::SoftKiller)

Pretty-print the SoftKiller grid configuration.
"""
function show(io::IO, sk::SoftKiller)
    if sk._ntotal <= 0
        print(io, "Uninitialized rectangular grid")
    else
        print(io,
              "rectangular grid with rapidity extent ",
              sk._ymin, " < rap < ", sk._ymax, "; ",
              "total tiles = ", sk._ntotal, "; ",
              "tile size Δy × Δφ = ", sk._dy, " × ", sk._dphi, ")")
    end
end

"""
    select_ABS_RAP_max(event, absrapmax)

Filter a collection of `PseudoJet`s, returning only those with absolute rapidity less than or equal to `absrapmax`.
"""
function select_ABS_RAP_max(event, absrapmax)
    filtered_events = filter(e -> begin
                                 abs(rapidity(e)) <= absrapmax
                             end, event)
    return filtered_events
end

"""
    softkiller(sk::SoftKiller, event::Vector{PseudoJet})

Apply the SoftKiller algorithm to an event (a vector of `PseudoJet`s).
Returns a tuple `(reduced_event, pt_threshold)`, where `reduced_event` is the filtered
event and `pt_threshold` is the computed pt threshold.
"""
function softkiller(sk::SoftKiller, event::Vector{PseudoJet})
    if (sk._ntotal < 2)
        throw("SoftKiller not properly initialized.")
    end

    # fills the vector of length n_tiles with 0's
    max_pt2 = fill(0.0, sk._ntotal)

    for ev in event
        index = tile_index(sk, ev)
        if (index < 0)
            continue
        end
        max_pt2[index] = max(max_pt2[index], pt2(ev))
    end

    pt2cut = median(max_pt2)

    reduced_event = filter(jet -> pt2(jet) >= pt2cut, event)

    pt_threshold = sqrt(pt2cut)

    return reduced_event, pt_threshold
end
