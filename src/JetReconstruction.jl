"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

# particle type definition
include("Particle.jl")
export energy, px, py, pz, pt, phi, mass, eta, kt, ϕ, η

# Philipp's pseudojet
include("pseudojet.jl")
export PseudoJet, rap, phi, pt2

# algorithmic part
include("Algo.jl")
export sequential_jet_reconstruct, kt_algo, anti_kt_algo, anti_kt_algo_alt, cambridge_aachen_algo

# jet serialisation (saving to file)
include("Serialize.jl")
export savejets, loadjets!, loadjets

# jet visualisation
include("JetVis.jl")
export jetsplot

end
