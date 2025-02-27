using PackageCompiler
import ArgParse
import JetReconstruction

function parse_args(args)
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--source-dir"
        help = "Directory containing the source files"
        arg_type = String
        default = splitdir(@__DIR__) |> first

        "--output-dir", "-o"
        help = "Directory to save the compiled library"
        arg_type = String
        default = joinpath(splitdir(@__DIR__) |> first, "JetReconstructionCompiled")
    end

    return ArgParse.parse_args(args, s)
end

function configure_file(template_path::String, output_path::String,
                        replacements::Dict{String, String})
    template = read(template_path, String)
    for (key, value) in replacements
        template = replace(template, "@$(key)@" => value)
    end
    open(output_path, "w") do io
        write(io, template)
    end
end

function (@main)(args)
    parsed_args = parse_args(args)
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

    cmake_input = joinpath(@__DIR__, "cmake", "JetReconstruction")
    cmake_output = joinpath(output_dir, "lib", "cmake", "JetReconstruction")

    version = pkgversion(JetReconstruction)
    cmake_project_version = "$(version.major).$(version.minor).$(version.patch)"

    @info "Copying CMake file to $cmake_output"
    mkpath(cmake_output)
    cp_file(input_dir, basename, output_dir) = cp(joinpath(input_dir, basename),
                                                  joinpath(output_dir, basename))
    cp_file(cmake_input, "JetReconstructionConfig.cmake", cmake_output)
    cp_file(cmake_input, "JetReconstructionTargets.cmake", cmake_output)
    configure_file(joinpath(cmake_input, "JetReconstructionConfigVersion.cmake.in"),
                   joinpath(cmake_output, "JetReconstructionConfigVersion.cmake"),
                   Dict("PROJECT_VERSION" => cmake_project_version))
end
