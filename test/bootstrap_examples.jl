# Bootstrap script to setup this local version of the JetReconstruction package
# for running the examples

println("Starting example boostrap script")

using Pkg
Pkg.instantiate()
Pkg.resolve()
Pkg.develop(path = "../JetReconstruction")
Pkg.status()

println("Finished example boostrap script")
