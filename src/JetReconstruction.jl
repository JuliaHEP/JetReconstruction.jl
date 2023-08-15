"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

# particle type definition
include("Particle.jl")
export energy, px, py, pz, pt, phi, mass, eta, kt, ϕ, η

# Philipp's pseudojet
include("Pseudojet.jl")
export PseudoJet, rap, phi, pt2

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

# Algorithmic part, tiled reconstruction strategy with SoA per tile
include("TiledAlgoSoATile.jl")
export tiled_jet_reconstruct_soa_tile

# Algorithmic part, tiled reconstruction strategy with global SoA
include("TiledAlgoSoAGlobal.jl")
export tiled_jet_reconstruct_soa_global

# Algorithmic part, tiled reconstruction strategy with linked list jet objects
include("TiledAlgoLL.jl")
export tiled_jet_reconstruct_ll

# jet serialisation (saving to file)
include("Serialize.jl")
export savejets, loadjets!, loadjets

# utility functions, useful for different primary scripts
include("Utils.jl")
export read_final_state_particles, pseudojets2vectors, final_jets

# jet visualisation
include("JetVis.jl")
export jetsplot

# JSON results
include("JSONresults.jl")
export FinalJet, FinalJets, JSON3

# Strategy to be used
@enum JetRecoStrategy Best N2Plain N2TiledLL N2TiledSoAGlobal N2TiledSoATile
export JetRecoStrategy, Best, N2Plain, N2TiledLL, N2TiledSoAGlobal, N2TiledSoATile

end
