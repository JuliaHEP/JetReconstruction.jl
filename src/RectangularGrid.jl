
using JetReconstruction

mutable struct RectangularGrid
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

    function RectangularGrid(rapmax::Float64, grid_size::Float64)
        grid = new(rapmax, -rapmax, grid_size, grid_size, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
        _setup_grid(grid)
        print(description(grid))

        return grid
    end

    function RectangularGrid(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)
        grid = new(rapmax, rapmin, drap, dphi, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
        _setup_grid(grid)
        print(description(grid))


        return grid
    end

    #Dummy constructor - gives unusable grid 
    RectangularGrid() = begin
        grid = RectangularGrid(-1.0, 1.0,
         -1.0, -1.0, -1, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
        _setup_grid(grid) #gives assertion error since invalid grid
        return grid 
    end

end

title_index(grid:: RectangularGrid, p::PseudoJet)::Int = begin 
    y_minus_ymin = rapidity(p) - grid._ymin
    if y_minus_ymin < 0 
        return -1   
    end 
    iy = y_minus_ymin * grid._inverse_dy
    if iy >= grid._nu
        return -1 
    end 

    iphi = phi(p)*grid._inverse_dphi
    if iphi == grid._nphi
        iphi =0
    end 

    iy*grid._nphi + iphi
end 

_setup_grid(grid:: RectangularGrid) = begin 
    @assert grid._ymax > grid._ymin
    @assert grid._requested_drap > 0
    @assert grid._requested_dphi > 0 

    ny_double = (grid._ymax-grid._ymin) / grid._requested_drap
    grid._ny = max(round(Int64, ny_double+0.5),1)
    grid._dy = (grid._ymax-grid._ymin) / grid._ny
    grid._inverse_dy = grid._ny/(grid._ymax-grid._ymin)

    grid._nphi = round(Int64,2* π / grid._requested_dphi + 0.5)
    grid._dphi = 2* π / grid._nphi
    grid._inverse_dphi = grid._nphi/2* π

    @assert grid._ny >=1 and grid._nphi >=1 

    grid._ntotal = grid._nphi * grid._ny
    grid._cell_area = grid._dy * grid._dphi

    #Selector implementation desn't exist 
end 
 
description(grid::RectangularGrid)::String = begin 
    #from definiton of is_initialised  in RectangularGrid.hh
    if grid._ntotal <= 0
        return "Uninitialised rectangular grid" 
    end 

    descr = "rectangular grid with rapidity extent $(grid._ymin) < rap < $(grid._ymax) \n total tiles  $(grid._ntotal) \n "
    descr *= "tile size drap x dphi = $(grid._dy) x $(grid._dphi)"

    #Selector implementation desn't exist 
    descr 
end 