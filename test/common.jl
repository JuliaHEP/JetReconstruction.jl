# This is a wrapper file that checks that we did not
# already include _common.jl.
#
# This is needed to allow sub-test scripts to be run
# independently, but to also allow the main test script
# to avoid double-including the common file.

if !@isdefined JETRECO_TEST_COMMON
    include("_common.jl")
end
