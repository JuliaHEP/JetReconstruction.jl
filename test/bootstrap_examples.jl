# Bootstrap script to setup this local version of the JetReconstruction package
# for running the examples

using Pkg
Pkg.instantiate()
Pkg.develop(path = "../JetReconstruction")
