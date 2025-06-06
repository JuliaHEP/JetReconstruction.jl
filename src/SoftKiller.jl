# Inspired by the SoftKiller implementation in the FastJet contrib package
# Original C++ code: https://fastjet.hepforge.org/contrib/

using JetReconstruction

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

    function SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)
        grid = new(rapmax, rapmin, drap, dphi, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
        _setup_grid!(grid)
        grid
    end

    function SoftKiller(rapmax::Float64, grid_size::Float64)
        grid = new(rapmax, -rapmax, grid_size, grid_size, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
                   0)
        _setup_grid!(grid)
        grid
    end
end

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

function select_ABS_RAP_max(event, absrapmax)
    filtered_events = filter(e -> begin
                                 abs(rapidity(e)) <= absrapmax
                             end, event)
    return filtered_events
end

function softkiller!(sk::SoftKiller, event::Vector{PseudoJet})
    if (sk._ntotal < 2)
        throw("SoftKiller not properly initialized.")
    end

    # fills the vector of length n_tiles with 0's
    max_pt2 = fill(0.0, sk._ntotal)

    for ev in event
        if (ev == isnothing)
            continue
        end
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
