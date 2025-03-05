using JetReconstruction

#do I need to account for different versions of fastjet? I assumed no 
mutable struct SoftKiller
    rapmin::Float64
    rapmax::Float64
    drap::Float64
    dphi::Float64
    grid::RectangularGrid
    #selector

end

#dummy SoftKiller - invalid 
function SoftKiller()
    RectangularGrid()
end

function SoftKiller(rapmin::Float64, rapmax::Float64, drap::Float64, dphi::Float64)
    print("4 variables \n")
    RectangularGrid(rapmin, rapmax, drap, dphi)

end

function SoftKiller(rapmax::Float64, grid_size::Float64)
    print("2 variables\n")
    RectangularGrid(rapmax, grid_size)
end

function set_up()

end 

function main()

end 

main()