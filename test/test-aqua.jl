# Use the Julia Aqua tests
# Aqua.jl: Auto QUality Assurance for Julia packages
using Aqua

include("common.jl")
@testset "Aqua.jl" begin
    Aqua.test_all(JetReconstruction)
end
