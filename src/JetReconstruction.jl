"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

# particle type definition
# include("Particle.jl")
# using .Particle
# export ParticleVector, PVec, energy, px, py, pz, pt, phi, mass, eta

# algorithmic part
include("Algo.jl")
using .Algo
export anti_kt, anti_kt!

end
