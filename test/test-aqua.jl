# Use the Julia Aqua tests
# Aqua.jl: Auto QUality Assurance for Julia packages
using Aqua

include("common.jl")
Aqua.test_all(JetReconstruction)
