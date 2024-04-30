"""
Jet reconstruction (reclustering) in Julia.
"""
module JetReconstruction

using LorentzVectorHEP

# Import from LorentzVectorHEP methods for those 4-vector types
pt2(p::LorentzVector) = LorentzVectorHEP.pt2(p)
phi(p::LorentzVector) = LorentzVectorHEP.phi(p)
rapidity(p::LorentzVector) = LorentzVectorHEP.rapidity(p)
px(p::LorentzVector) = LorentzVectorHEP.px(p)
py(p::LorentzVector) = LorentzVectorHEP.py(p)
pz(p::LorentzVector) = LorentzVectorHEP.pz(p)
energy(p::LorentzVector) = LorentzVectorHEP.energy(p)

pt2(p::LorentzVectorCyl) = LorentzVectorHEP.pt2(p)
phi(p::LorentzVectorCyl) = LorentzVectorHEP.phi(p)
rapidity(p::LorentzVectorCyl) = LorentzVectorHEP.rapidity(p)
px(p::LorentzVectorCyl) = LorentzVectorHEP.px(p)
py(p::LorentzVectorCyl) = LorentzVectorHEP.py(p)
pz(p::LorentzVectorCyl) = LorentzVectorHEP.pz(p)
energy(p::LorentzVectorCyl) = LorentzVectorHEP.energy(p)

# Pseudojet type
include("Pseudojet.jl")
export PseudoJet

# Jet reconstruction strategies and algorithms (enums!)
include("AlgorithmStrategyEnums.jl")
export JetRecoStrategy, JetAlgorithm

# ClusterSequence type
include("ClusterSequence.jl")
export ClusterSequence, inclusive_jets, exclusive_jets

## N2Plain algorithm
# Algorithmic part for simple sequential implementation
include("PlainAlgo.jl")
export plain_jet_reconstruct

## N2Tiled algorithm
# Common pieces
include("TiledAlgoUtils.jl")
# Algorithmic part, tiled reconstruction strategy with linked list jet objects
include("TiledAlgoLL.jl")
export tiled_jet_reconstruct

## Generic algorithm, which can switch strategy dynamically
include("GenericAlgo.jl")
export jet_reconstruct

# Simple HepMC3 reader
include("HepMC3.jl")

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

end
