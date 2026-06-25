using PackageCompiler
import ArgParse
import JetReconstruction
import JuliaC

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

        "--juliac"
        help = "Use JuliaC.jl for compilation instead of PackageCompiler.jl"
        action = :store_true
    end

    return ArgParse.parse_args(args, s)
end

function cp_file(input_dir, basename, output_dir)
    cp(joinpath(input_dir, basename),
       joinpath(output_dir, basename); force = true)
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

const C_ENUM_SPECS = ((placeholder = "JETRECONSTRUCTION_JETALGORITHM_ENUM",
                       enum_type = JetReconstruction.JetAlgorithm.Algorithm,
                       prefix = "JETRECONSTRUCTION_JETALGORITHM"),
                      (placeholder = "JETRECONSTRUCTION_RECOSTRATEGY_ENUM",
                       enum_type = JetReconstruction.RecoStrategy.Strategy,
                       prefix = "JETRECONSTRUCTION_RECOSTRATEGY"),
                      (placeholder = "JETRECONSTRUCTION_RECOMBINATIONSCHEME_ENUM",
                       enum_type = JetReconstruction.RecombinationScheme.Recombine,
                       prefix = "JETRECONSTRUCTION_RECOMBINATIONSCHEME"))

function render_c_enum_entries(enum_type::Type{<:Enum}, prefix::AbstractString)
    vals = collect(instances(enum_type))
    join(("  $(prefix)_$(uppercase(String(Symbol(v)))) = $(Int(v))" *
          (i == lastindex(vals) ? "" : ",")
          for (i, v) in pairs(vals)),
         "\n")
end

function c_enum_replacements()
    return Dict(spec.placeholder => render_c_enum_entries(spec.enum_type, spec.prefix)
                for spec in C_ENUM_SPECS)
end

function compile_w_packagecompiler(source_dir, output_dir)
    return @elapsed PackageCompiler.create_library(source_dir, output_dir;
                                                   lib_name = "jetreconstruction",
                                                   precompile_execution_file = [joinpath(@__DIR__,
                                                                                         "precompile_execution.jl")],
                                                   incremental = false,
                                                   filter_stdlibs = true,
                                                   force = true)
end

function compile_w_juliac(source_dir, output_dir)
    img = JuliaC.ImageRecipe(output_type = "--output-lib",
                             trim_mode = "no",
                             file = joinpath(source_dir, "src", "JetReconstruction.jl"),
                             project = source_dir,
                             add_ccallables = true,
                             verbose = false)

    link = JuliaC.LinkRecipe(image_recipe = img,
                             outname = joinpath(output_dir, "lib", "libjetreconstruction"))

    bun = JuliaC.BundleRecipe(link_recipe = link,
                              output_dir = nothing)
    return @elapsed begin
        JuliaC.compile_products(img)
        JuliaC.link_products(link)
        JuliaC.bundle_products(bun)
    end
end

function (@main)(args)
    parsed_args = parse_args(args)
    source_dir = parsed_args["source-dir"]
    output_dir = parsed_args["output-dir"]

    @info "Compiling package from $source_dir"
    @info "Creating library in $output_dir"

    compilation_time = if parsed_args["juliac"]
        @info "Using JuliaC.jl"
        compile_w_juliac(source_dir, output_dir)
    else
        @info "Using PackageCompiler.jl"
        compile_w_packagecompiler(source_dir, output_dir)
    end
    @info "Compiled in $(compilation_time) seconds"
    compiler = parsed_args["juliac"] ? "JETRECONSTRUCTION_COMPILER_JULIAC" :
               "JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER"

    includes_input = joinpath(@__DIR__, "include")
    includes_output = joinpath(output_dir, "include")
    mkpath(includes_output)
    @info "Copying header files to $includes_output"
    configure_file(joinpath(includes_input, "JetReconstruction.h.in"),
                   joinpath(includes_output, "JetReconstruction.h"),
                   merge(Dict("JETRECONSTRUCTION_COMPILER" => compiler),
                         c_enum_replacements()))

    cmake_input = joinpath(@__DIR__, "cmake", "JetReconstruction")
    cmake_output = joinpath(output_dir, "lib", "cmake", "JetReconstruction")
    @info "Copying CMake files to $cmake_output"
    mkpath(cmake_output)

    version = pkgversion(JetReconstruction)
    cmake_project_version = "$(version.major).$(version.minor).$(version.patch)"

    configure_file(joinpath(cmake_input, "JetReconstructionConfig.cmake.in"),
                   joinpath(cmake_output, "JetReconstructionConfig.cmake"),
                   Dict("JETRECONSTRUCTION_COMPILER" => compiler))
    configure_file(joinpath(cmake_input, "JetReconstructionConfigVersion.cmake.in"),
                   joinpath(cmake_output, "JetReconstructionConfigVersion.cmake"),
                   Dict("PROJECT_VERSION" => cmake_project_version))
    cp_file(cmake_input, "JetReconstructionTargets.cmake", cmake_output)
    return 0
end
