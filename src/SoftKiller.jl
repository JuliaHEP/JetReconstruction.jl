using Statistics

"""
    struct SoftKiller{T <: Real}

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
- `_ymin::T`: Minimum rapidity of the grid.
- `_ymax::T`: Maximum rapidity of the grid.
- `_requested_drap::T`: Requested grid spacing in rapidity.
- `_requested_dphi::T`: Requested grid spacing in phi.
- `_ntotal::Int`: Total number of tiles.
- `_dy::T`: Actual grid spacing in rapidity.
- `_dphi::T`: Actual grid spacing in phi.
- `_cell_area::T`: Area of a single tile.
- `_inverse_dy::T`: Inverse of rapidity grid spacing.
- `_inverse_dphi::T`: Inverse of phi grid spacing.
- `_ny::Int`: Number of tiles in rapidity.
- `_nphi::Int`: Number of tiles in phi.

# Constructors
- `SoftKiller(rapmin::T, rapmax::T, drap::T, dphi::T)`: 
  Construct a grid from `rapmin` to `rapmax` in rapidity, with tile sizes `drap` and `dphi`.
- `SoftKiller(rapmax::T, grid_size::T)`: 
  Construct a square grid from `-rapmax` to `rapmax` in rapidity, with tile size `grid_size`.

"""
struct SoftKiller{T <: Real}
    _ymin::T
    _ymax::T
    _requested_drap::T
    _requested_dphi::T
    _ntotal::Int
    _dy::T
    _dphi::T
    _cell_area::T
    _inverse_dy::T
    _inverse_dphi::T
    _ny::Int
    _nphi::Int

    """
        SoftKiller(rapmin::T, rapmax::T, drap::T, dphi::T)

    Construct a SoftKiller grid from `rapmin` to `rapmax` in rapidity, with tile sizes `drap` and `dphi`.
    """

    function SoftKiller(rapmin::T, rapmax::T, drap::T, dphi::T) where {T <: Real}
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

        new{T}(rapmin, rapmax, drap, dphi, ntotal, dy, dphi_final, cell_area,
               inverse_dy, inverse_dphi, ny, nphi)
    end

    """
        SoftKiller(rapmax::T, grid_size::T)

    Construct a square SoftKiller grid from `-rapmax` to `rapmax` in rapidity, with tile size `grid_size`.
    """
    function SoftKiller(rapmax::T, grid_size::T) where {T <: Real}
        return SoftKiller(-rapmax, rapmax, grid_size, grid_size)
    end
end

"""
    tile_index(sk::SoftKiller{U}, p::PseudoJet{T}) where {U, T}

Return the tile index for a given `PseudoJet` `p` in the SoftKiller grid `sk`.
Returns -1 if the jet is outside the grid bounds.
"""
function tile_index(sk::SoftKiller{U}, p::PseudoJet{T}) where {U <: Real, T <: Real}
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
    show(io::IO, sk::SoftKiller{U}) where {U}

Pretty-print the SoftKiller grid configuration.
"""
function show(io::IO, sk::SoftKiller{U}) where {U <: Real}
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
    softkiller(sk::SoftKiller{U}, event::Vector{J}) where {U, T <: Real, J <: PseudoJet{T}}

Apply the SoftKiller algorithm to an event (a vector of `PseudoJet`s).
Returns a tuple `(reduced_event, pt_threshold)`, where `reduced_event` is the filtered
event and `pt_threshold` is the computed pt threshold.
"""
function softkiller(sk::SoftKiller{U},
                    event::Vector{J}) where {U, T <: Real, J <: PseudoJet{T}}
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

    reduced_event = Vector{J}()
    sizehint!(reduced_event, length(event))
    hist_index = 1
    for jet in event
        if pt2(jet) >= pt2cut
            push!(reduced_event, PseudoJet{T}(jet; cluster_hist_index = hist_index))
            hist_index += 1
        end
    end

    pt_threshold = sqrt(pt2cut)

    return reduced_event, pt_threshold
end
