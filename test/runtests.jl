#! /usr/bin/env julia

include("common.jl")

function main()
    # Algorithm/power consistency checks
    include("test-algpower-consistency.jl")

    # Jet utilities tests
    include("test-jet-utils.jl")

    # jet_reconstruct() interface check
    include("test-jet_reconstruct-interface.jl")

    # New test structure, factorised tests for pp and e+e-
    include("test-pp-reconstruction.jl")
    include("test-ee-reconstruction.jl")

    # Test jets selection
    include("test-selection.jl")

    # Compare inputting data in PseudoJet with using a LorentzVector
    include("test-compare-types.jl")

    # Check jet constituents
    include("test-constituents.jl")

    # Suppress these tests for now, as the examples Project.toml is rather heavy
    # because of the GLMakie dependency, plus on a CI there is no GL subsystem,
    # so things fail. The examples should be restructured to have a cleaner set
    # of examples that are good to run in the CI.
    #
    # Now run a few tests with our examples
    # include("tests_examples.jl")

    # Utility support tests
    include("test-utils.jl")

    # Substructure tests
    include("test-substructure.jl")

    # Test with Aqua (https://juliatesting.github.io/Aqua.jl/stable/)
    include("test-aqua.jl")
end

logger = ConsoleLogger(stdout, Logging.Warn)
global_logger(logger)
main()
