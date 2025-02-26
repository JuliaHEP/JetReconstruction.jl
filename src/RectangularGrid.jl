#RectangularGrid 
#using Pkg
#using Pkg
using JetReconstruction
#Pkg.add("LorentzVectorHEP")
#kg.add("MuladdMacro")
#Pkg.add("StructArrays")

#Pkg.add("PseudoJet.jl")
#include("PseudoJet.jl") 
#Pkg.develop(path="/Users/emadimtrova/Desktop/uni/spring25/UTRA/JetReconstruction.jl/src/PseudoJet.jl")
#Pkg.add("PseudoJet")
#using Pseudojet
#include("fastjet/Pseudojet.jl")

#using Pkg
#Pkg.develop(path="/Users/emadimtrova/Desktop/uni/spring25/UTRA/JetReconstruction.jl/src/")
#include("Pseudojet.jl")  # This includes JetReconstruct.jl
#using Pseudojet

mutable struct RectangularGrid
    _ymax::Float64
    _ymin::Float64
    _requested_drap::Float64
    _requested_dphi::Float64
    _ntotal::Int
    _ngood::Int
    _dy::Float64
    _dphi::Float64
    _cell_area::Float64
    _inverse_dy::Float64
    _inverse_dphi::Float64
    _ny::Int
    _nphi::Int
end

#Dummy constructor 
function RectangularGrid()
    RectangularGrid(-1.0, 1.0, -1.0, -1.0, -1,   
        -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)
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
    grid._ny = max(int(ny_double+0.5),1)
    grid._dy = (grid._ymax-grid._ymin) / grid._ny
    grid._inverse_dy = grid._ny/(grid._ymax-grid._ymin)

    grid._nphi = 2* π / grid._requested_dphi + 0.5
    grid._dphi = 2* π / _nphi
    grid._inverse_dphi = grid._nphi/2* π

    @assert grid._ny >=1 and grid._nphi >=1 

    grid._ntotal = grid._nphi * grid._ny
    grid._cell_area = grid._dy * grid._dphi

    #Selector implementation desn't exist 
end 
 
description(grid::RectangularGrid)::string = begin 
    #from definiton of is_initialised  in RectangularGrid.hh
    if grid._ntotal <= 0
        return "Uninitialised rectangular grid" 
    end 

    descr = "rectangular grid with rapidity extent $(grid._ymin) < rap < $(grid._ymax), "
    descr *= "tile size drap x dphi = $(grid._dy) x $(grid._dphi)"

    #Selector implementation desn't exist 

    descr 
end 

