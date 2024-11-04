#! /usr/bin/env julia

using Test

const example_dir = joinpath(pkgdir(JetReconstruction), "examples")
const test_dir = joinpath(pkgdir(JetReconstruction), "test")

# Need to ensure examples is properly instantiated
run(`julia --project=$example_dir $test_dir/bootstrap_examples.jl`)

# Now run a few tests with our examples
@testset "Jet Reconstruction Examples" begin
    @test success(run(pipeline(`julia --project=$example_dir $example_dir/jetreco.jl -S N2Plain -A AntiKt -R 1.0 $events_file`,
                               devnull)))
    @test success(run(pipeline(`julia --project=$example_dir $example_dir/jetreco.jl -S N2Tiled -p 1 -R 1.0 $events_file`,
                               devnull)))
    @test success(run(pipeline(`julia --project=$example_dir $example_dir/jetreco.jl -p 1.5 -A GenKt -R 1.0 $events_file`,
                               devnull)))

    @test success(run(pipeline(`julia --project=$example_dir $example_dir/instrumented-jetreco.jl -S N2Plain -A AntiKt -R 1.0 $events_file`,
                               devnull)))
    @test success(run(pipeline(`julia --project=$example_dir $example_dir/instrumented-jetreco.jl -S N2Tiled -p 1 -R 1.0 $events_file`,
                               devnull)))
    @test success(run(pipeline(`julia --project=$example_dir $example_dir/instrumented-jetreco.jl -p 1.5 -A GenKt -R 1.0 $events_file`,
                               devnull)))

    @test success(run(pipeline(`julia --project=$example_dir $example_dir/jetreco-constituents.jl`,
                               devnull)))
    @test success(run(pipeline(`julia --project=$example_dir $example_dir/visualise-jets.jl -A AntiKt --event 5 -R 2.0 $events_file $example_dir/jetvis.png`,
                               devnull)))
end
