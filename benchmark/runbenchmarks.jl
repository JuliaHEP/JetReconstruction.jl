"""Benchmark runner for JetReconstruction.jl package."""

using ArgParse
using BenchmarkTools
using PkgBenchmark

"""
    find_trials(group::BenchmarkGroup; path=String[])

Recursively traverse a BenchmarkGroup and return a flat dictionary of all BenchmarkTools.Trial objects found
with their hierarchical path as the key.
"""
function find_trials(group::BenchmarkGroup; path = String[])
    leaves = Dict{String, BenchmarkTools.Trial}()

    for (k, v) in group
        newpath = [path..., string(k)]
        if v isa BenchmarkGroup
            merge!(leaves, find_trials(v; path = newpath))
        elseif v isa BenchmarkTools.Trial
            key = join(newpath, "/")
            leaves[key] = v
        end
    end

    return leaves
end

function display_results(flat_results::Dict)
    for key in sort(collect(keys(flat_results)))
        println("• $key")
        display(flat_results[key])
    end
end

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true,
                         description = "Run JetReconstruction.jl benchmark suite")
    @add_arg_table! s begin
        "--verbose", "-v"
        help = "Enable verbose output, print results to console"
        action = :store_true

        "--threads", "-t"
        help = "Number of threads for the benchmark process (forwarded via BenchmarkConfig)"
        arg_type = Int
        default = 1

        "output"
        help = "Write benchmark results to BenchmarkTools.BenchmarkGroup JSON formatted file"
    end
    return parse_args(args, s)
end

function (@main)(args)
    parsed_args = parse_command_line(args)
    n_threads = parsed_args["threads"]
    @info "Benchmark process threads: $(n_threads)"
    config = PkgBenchmark.BenchmarkConfig(juliaargs = ["-t", string(n_threads)])
    results = benchmarkpkg(dirname(@__DIR__); config = config,
                           resultfile = parsed_args["output"])
    if parsed_args["verbose"]
        (display_results ∘ find_trials ∘ PkgBenchmark.benchmarkgroup)(results)
    end
    return 0
end

@isdefined(var"@main") ? (@main) : exit(main(ARGS))
