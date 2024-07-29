# Bootstrap script to setup this local version of the JetReconstruction package
# for running the examples

println("Starting example boostrap script")

if isfile(joinpath(@__DIR__, "..", "examples", "Manifest.toml"))
    println("Exisiting Manifest.toml found - assuming environment is already setup")
    exit(0)
end

if isdir(joinpath("..", "JetReconstruction.jl"))
    symlink("JetReconstruction.jl", "../JetReconstruction")
else

println("Trying local path: $local_path (from $(pwd()))")

using Pkg
Pkg.instantiate()
Pkg.resolve()
Pkg.develop(path = "../JetReconstruction")
Pkg.status()

println("Finished example boostrap script")
