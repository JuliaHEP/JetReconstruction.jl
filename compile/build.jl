using PackageCompiler
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--source-dir"
    help = "Directory containing the source files"
    arg_type = String
    default = splitdir(@__DIR__) |> first

    "--output-dir"
    help = "Directory to save the compiled library"
    arg_type = String
    default = joinpath(splitdir(@__DIR__) |> first, "JetReconstructionCompiled")
end

parsed_args = parse_args(s)
source_dir = parsed_args["source-dir"]
output_dir = parsed_args["output-dir"]

@info "Compiling package from $source_dir"
@info "Creating library in $output_dir"
PackageCompiler.create_library(source_dir, output_dir;
                               lib_name = "jetreconstruction",
                               header_files = [joinpath(@__DIR__, "include",
                                                        "JetReconstruction.h")],
                               precompile_execution_file = [joinpath(@__DIR__,
                                                                     "precompile_execution.jl")],
                               incremental = false,
                               filter_stdlibs = true,
                               force = true)
