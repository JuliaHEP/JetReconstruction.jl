# Use the Julia Aqua tests
# Aqua.jl: Auto QUality Assurance for Julia packages
include("common.jl")

using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(JetReconstruction)
end
