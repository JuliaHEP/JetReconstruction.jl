#! /usr/bin/env julia

include("common.jl")

@testset "Basic jetreco.jl reconstruction examples" begin
    @test success(run(pipeline(`julia --project=examples $(@__DIR__)/../examples/jetreco.jl --algorithm=CA -R 0.4 $(@__DIR__)/../test/data/events.pp13TeV.hepmc3.zst`,
                               devnull)))
    @test success(run(pipeline(`julia --project=examples $(@__DIR__)/../examples/jetreco.jl --algorithm=EEKt -R 1.0 -p -1.0 $(@__DIR__)/../test/data/events.eeH.hepmc3.zst`,
                               devnull)))
    @test success(run(pipeline(`julia --project=examples $(@__DIR__)/../examples/jetreco.jl --algorithm=Durham $(@__DIR__)/../test/data/events.eeH.hepmc3.zst --dump test-output-ee.csv`,
                               devnull)))
    rm("test-output-ee.csv"; force = true)
end
