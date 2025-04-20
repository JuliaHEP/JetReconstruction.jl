
using JetReconstruction

#abstract type RectangularGrid end

mutable struct RectangularGrid <: TilingBase
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
        grid = new(rapmax, -rapmax, grid_size, grid_size, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
                   0)
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
    #RectangularGrid() = begin
    #    grid = RectangularGrid(-1.0, 1.0,
    #     -1.0, -1.0, -1, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
    #    _setup_grid(grid) #gives assertion error since invalid grid
    #   return grid 
    #end

end
