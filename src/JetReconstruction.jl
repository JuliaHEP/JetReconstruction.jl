"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

using LorentzVectorHEP

# particle type definition
# include("Particle.jl")
# export energy, px, py, pz, pt, phi, mass, eta, kt, ϕ, η

# Philipp's pseudojet
include("Pseudojet.jl")
## As this is an internal EDM class, we perhaps shouldn't export this stuff...
# export PseudoJet, rap, phi, pt, pt2, px, py, pz, pt, phi, mass, eta

# Simple HepMC3 reader
include("HepMC3.jl")
export HepMC3

## N2Plain algorithm
# Algorithmic part for simple sequential implementation
include("Algo.jl")
export sequential_jet_reconstruct, kt_algo, anti_kt_algo, anti_kt_algo_alt, cambridge_aachen_algo

## Tiled algorithms
# Common pieces
include("TiledAlgoUtils.jl")

# Algorithmic part, tiled reconstruction strategy with linked list jet objects
include("TiledAlgoLL.jl")
export tiled_jet_reconstruct_ll

# jet serialisation (saving to file)
include("Serialize.jl")
export savejets, loadjets!, loadjets

# utility functions, useful for different primary scripts
include("Utils.jl")
export read_final_state_particles, read_final_state_particles_lv, pseudojets2vectors, final_jets

# jet visualisation
include("JetVis.jl")
export jetsplot

# JSON results
include("JSONresults.jl")
export FinalJet, FinalJets, JSON3

# Strategy to be used
## Maybe an enum is not the best idea, use type dispatch instead?
@enum JetRecoStrategy Best N2Plain N2Tiled
export JetRecoStrategy, Best, N2Plain, N2Tiled

end
