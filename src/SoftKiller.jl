using JetReconstruction

"""
SoftKiller

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
- `_ntotal::Int64`: Total number of tiles.
- `_ngood::Int64`: Number of tiles with at least one particle.
- `_dy::Float64`: Actual grid spacing in rapidity.
- `_dphi::Float64`: Actual grid spacing in phi.
- `_cell_area::Float64`: Area of a single tile.
- `_inverse_dy::Float64`: Inverse of rapidity grid spacing.
- `_inverse_dphi::Float64`: Inverse of phi grid spacing.
- `_ny::Int64`: Number of tiles in rapidity.
- `_nphi::Int64`: Number of tiles in phi.

# Constructors
- `SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)`: 
  Construct a grid from `rapmin` to `rapmax` in rapidity, with tile sizes `drap` and `dphi`.
- `SoftKiller(rapmax::Float64, grid_size::Float64)`: 
  Construct a square grid from `-rapmax` to `rapmax` in rapidity, with tile size `grid_size`.

"""
mutable struct SoftKiller
    _ymax::Float64
    _ymin::Float64
    _requested_drap::Float64
    _requested_dphi::Float64
    _ntotal::Int64
    _ngood::Int64
    _dy::Float64
    _dphi::Float64
    _cell_area::Float64
    _inverse_dy::Float64
    _inverse_dphi::Float64
    _ny::Int64
    _nphi::Int64

    """
        SoftKiller(rapmax::Float64, grid_size::Float64)

    Construct a square SoftKiller grid from `-rapmax` to `rapmax` in rapidity, with tile size `grid_size`.
    """
    function SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)
        grid = new(rapmax, rapmin, drap, dphi, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
        _setup_grid!(grid)
        grid
    end

    """
        SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)

    Construct a SoftKiller grid from `rapmin` to `rapmax` in rapidity, with tile sizes `drap` and `dphi`.
    """
    function SoftKiller(rapmax::Float64, grid_size::Float64)
        grid = new(rapmax, -rapmax, grid_size, grid_size, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
                   0)
        _setup_grid!(grid)
        grid
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

    iy = round(Int64, y_minus_ymin * sk._inverse_dy)
    if iy >= sk._ny
        return -1
    end

    iphi = round(Int64, phi(p) * sk._inverse_dphi)
    if iphi == sk._nphi
        iphi = 0
    end

    res = round(Int64, iy * sk._nphi + iphi)

    res + 1
end

"""
    _setup_grid!(sk::SoftKiller)

Internal function to initialize the grid parameters for a `SoftKiller` instance.
"""
function _setup_grid!(sk::SoftKiller)
    @assert sk._ymax > sk._ymin
    @assert sk._requested_drap > 0
    @assert sk._requested_dphi > 0

    ny_double = (sk._ymax - sk._ymin) / sk._requested_drap
    sk._ny = max(round(Int64, ny_double + 0.5), 1)
    sk._dy = (sk._ymax - sk._ymin) / sk._ny
    sk._inverse_dy = sk._ny / (sk._ymax - sk._ymin)

    sk._nphi = round(Int64, (2 * π) / sk._requested_dphi + 0.5)
    sk._dphi = (2 * π) / sk._nphi
    sk._inverse_dphi = sk._nphi / (2 * π)

    @assert sk._ny>=1 and sk._nphi>=1

    sk._ntotal = sk._nphi * sk._ny
    sk._cell_area = sk._dy * sk._dphi
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
    softkiller!(sk::SoftKiller, event::Vector{PseudoJet})

Apply the SoftKiller algorithm to an event (a vector of `PseudoJet`s).
Returns a tuple `(reduced_event, pt_threshold)`, where `reduced_event` is the filtered
event and `pt_threshold` is the computed pt threshold.
"""
function softkiller!(sk::SoftKiller, event::Vector{PseudoJet})
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

    sort!(max_pt2)

    int_median_pos = length(max_pt2) ÷ 2
    pt2cut = (1 + 1e-12) * max_pt2[int_median_pos]

    indices = Int64[]
    for (i, ps_jet) in enumerate(event)
        if ps_jet === nothing || pt2(ps_jet) >= pt2cut
            push!(indices, i)
        end
    end

    reduced_event = event[indices]

    pt_threshold = sqrt(pt2cut)

    return reduced_event, pt_threshold
end
