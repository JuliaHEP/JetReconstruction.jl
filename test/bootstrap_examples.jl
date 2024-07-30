# Bootstrap script to setup this local version of the JetReconstruction package
# for running the examples

println("Starting example boostrap script")

if isfile(joinpath(@__DIR__, "..", "examples", "Manifest.toml"))
    println("Exisiting Manifest.toml found - assuming environment is already setup")
    exit(0)
end

local_pkg_path = realpath(joinpath(@__DIR__, "..", "..", "JetReconstruction"))

if !isdir(local_pkg_path)
    # Try a symlink to the current checkout
    local_checkout_path = realpath(joinpath(@__DIR__, "..", ".."))
    println("Creating symlink from $local_pkg_path to $local_checkout_path")
    symlink($local_checkout_path, local_pkg_path)
end

println("Setting up examples project")

using Pkg
Pkg.instantiate()
Pkg.resolve()
println(("Seting JetReconstruction development package path: $local_pkg_path"))
Pkg.develop(path = local_pkg_path)
Pkg.status()

println("Finished examples boostrap script")
