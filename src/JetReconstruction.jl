"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

# particle type definition
include("Particle.jl")
using .Particle
export energy, px, py, pz, pt, phi, mass, eta, kt, ϕ, η

# algorithmic part
include("Algo.jl")
using .Algo
export anti_kt, anti_kt_alt

end
