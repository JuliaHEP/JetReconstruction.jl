using PackageCompiler

output_directory = joinpath(splitdir(@__DIR__) |> first, "JetReconstructionCompiled")

@info "Creating library in $output_directory"
PackageCompiler.create_library(".", output_directory;
                               lib_name = "jetreconstruction",
                               header_files = [joinpath(@__DIR__, "include",
                                                        "JetReconstruction.h")],
                               #       precompile_execution_file = [joinpath(@__DIR__,
                               #                                             "precompile_execution.jl")],
                               #       precompile_statements_file = [jointpath(@__DIR__,
                               #                                               "precompile_statements.jl")],
                               incremental = false,
                               filter_stdlibs = true,
                               force = true)
