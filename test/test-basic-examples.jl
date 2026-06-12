#! /usr/bin/env julia

include("common.jl")

const JULIA_CMD = Base.julia_cmd()
const EXAMPLES_PROJECT = joinpath(@__DIR__, "..", "examples")

@testset "Instantiate examples environment" begin
    instantiate_examples = """
        using Pkg
        cd($(repr(EXAMPLES_PROJECT)))
        Pkg.activate($(repr(EXAMPLES_PROJECT)))
        Pkg.develop(path="..")
        Pkg.instantiate()
    """
    @test success(run(`$JULIA_CMD -e $instantiate_examples`))
end

@testset "Basic jetreco.jl reconstruction examples" begin
    @test success(run(pipeline(`$JULIA_CMD --project=$EXAMPLES_PROJECT $(@__DIR__)/../examples/jetreco.jl --algorithm=CA -R 0.4 $(@__DIR__)/../test/data/events.pp13TeV.hepmc3.zst`,
                               devnull)))
    @test success(run(pipeline(`$JULIA_CMD --project=$EXAMPLES_PROJECT $(@__DIR__)/../examples/jetreco.jl --algorithm=EEKt -R 1.0 -p -1.0 $(@__DIR__)/../test/data/events.eeH.hepmc3.zst`,
                               devnull)))
    @test success(run(pipeline(`$JULIA_CMD --project=$EXAMPLES_PROJECT $(@__DIR__)/../examples/jetreco.jl --algorithm=Durham $(@__DIR__)/../test/data/events.eeH.hepmc3.zst --dump test-output-ee.csv`,
                               devnull)))
    rm("test-output-ee.csv"; force = true)
end

@testset "Basic instrumented-jetreco.jl reconstruction example" begin
    @test success(run(pipeline(`$JULIA_CMD --project=$EXAMPLES_PROJECT $(@__DIR__)/../examples/instrumented-jetreco.jl --maxevents=1 --nsamples=16 --algorithm=CA -R 0.4 --plot $(@__DIR__)/../test/data/events.pp13TeV.hepmc3.zst`,
                               devnull)))
end
