"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

# particle type definition
include("Particle.jl")
export energy, px, py, pz, pt, phi, mass, eta, kt, ϕ, η

# algorithmic part
include("Algo.jl")
export kt_algo, anti_kt_algo, anti_kt_algo_alt

# jet serialisation (saving to file)
include("Serialize.jl")
export savejets, loadjets!, loadjets

# jet visualisation
include("JetVis.jl")
export jetsplot

end
