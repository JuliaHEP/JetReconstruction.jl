#! /usr/bin/env julia

include("common.jl")

function main()
    # Algorithm/power consistency checks
    include("test-algpower-consistency.jl")

    # Basic tests for the Jet types
    include("test-jet-types.jl")
    include("test-jet-constructors.jl")
    include("test-edm4hep.jl")

    # Jet utilities tests
    include("test-jet-utils.jl")

    # jet_reconstruct() interface check
    include("test-jet-reconstruct-interface.jl")

    # Core tests of EE and PP reconstruction, with comparison to Fastjet
    include("test-pp-reconstruction.jl")
    include("test-ee-reconstruction.jl")

    # Valencia algorithm tests
    include("test-valencia.jl")

    # Test jets selection
    include("test-selection.jl")

    # Compare inputting data in PseudoJet with using a LorentzVector
    include("test-compare-types.jl")

    # Check jet constituents
    include("test-constituents.jl")

    # Utility support tests
    include("test-utils.jl")

    # Check jet selectors to different jet types
    include("test-jet-selectors.jl")

    # Substructure tests
    include("test-substructure.jl")

    # LundPlane tests
    include("test-lundplane.jl")

    # SoftKiller tests
    include("test-softkiller.jl")

    # Test that reconstruction runs with different numerical precision
    include("test-reco-numtypes.jl")

    # Test basic examples (where the dependencies are not too heavy)
    include("test-basic-examples.jl")

    # Test with Aqua (https://juliatesting.github.io/Aqua.jl/stable/)
    include("test-aqua.jl")
end

logger = ConsoleLogger(stdout, Logging.Warn)
global_logger(logger)
main()
