# EnumX provides scoped enums, which are nicer
using EnumX

"""
    enum RecoStrategy

Scoped enumeration (using EnumX) representing the different strategies for jet reconstruction.

## Fields
- `Best`: The best strategy.
- `N2Plain`: The plain N2 strategy.
- `N2Tiled`: The tiled N2 strategy.
"""
@enumx T=Strategy RecoStrategy Best N2Plain N2Tiled
const AllJetRecoStrategies = [String(Symbol(x)) for x in instances(RecoStrategy.Strategy)]

"""
    enum T

Scoped enumeration (using EnumX) representing different jet algorithms used in the JetReconstruction module.

## Fields
- `AntiKt`: The Anti-Kt algorithm.
- `CA`: The Cambridge/Aachen algorithm.
- `Kt`: The Inclusive-Kt algorithm.
- `GenKt`: The Generalised Kt algorithm (with arbitrary power).
- `EEKt`: The Generalised e+e- kt algorithm.
- `Durham`: The e+e- kt algorithm, aka Durham.
- `Valencia`: The Valencia e+e- algorithm.
"""
@enumx T=Algorithm JetAlgorithm AntiKt CA Kt GenKt EEKt Durham Valencia
const AllJetRecoAlgorithms = [String(Symbol(x)) for x in instances(JetAlgorithm.Algorithm)]

"""
    const varpower_algorithms

A constant array that contains the jet algorithms for which power is variable.

"""
const varpower_algorithms = [JetAlgorithm.GenKt, JetAlgorithm.EEKt, JetAlgorithm.Valencia]

"""
    algorithm2power

A dictionary that maps algorithm names to their corresponding power values.
"""
const algorithm2power = Dict(JetAlgorithm.AntiKt => -1,
                             JetAlgorithm.CA => 0,
                             JetAlgorithm.Kt => 1,
                             JetAlgorithm.Durham => 1)

"""
    get_algorithm_power(; algorithm::JetAlgorithm.Algorithm, p::Union{Real, Nothing}) -> Real

Get the algorithm power

This function returns appropriate power value for the specified jet algorithm if the algorithm
is a fixed power algorithm (like `AntiKt`, `CA`, `Kt`, or `Durham`).
If the algorithm is generalized (like `GenKt` or `EEKt`), it requires a power value
to be specified and the function returns the same value.

# Arguments
- `p::Union{Real, Nothing}`: The power value.
- `algorithm::JetAlgorithm.Algorithm`: The algorithm.

# Returns
Resolved algorithm power value.

# Throws
- `ArgumentError`: If no power is specified for a generalized algorithm

"""
function get_algorithm_power(; algorithm::JetAlgorithm.Algorithm, p::Union{Real, Nothing})
    # For the case where an algorithm has a variable power value
    # we need to check that the power value and algorithm are both specified
    if algorithm in varpower_algorithms
        if isnothing(p)
            throw(ArgumentError("Power must be specified for algorithm $algorithm"))
        end
        # And then we just return what these are...
        return p
    end
    # Otherwise the algorithm has a fixed power value
    return algorithm2power[algorithm]
end

"""
    is_pp(algorithm::JetAlgorithm.Algorithm)

Check if the algorithm is a pp reconstruction algorithm.

# Returns
`true` if the algorithm is a pp reconstruction algorithm, `false` otherwise.
"""
function is_pp(algorithm::JetAlgorithm.Algorithm)
    return algorithm in [
        JetAlgorithm.AntiKt,
        JetAlgorithm.CA,
        JetAlgorithm.Kt,
        JetAlgorithm.GenKt
    ]
end

"""
    is_ee(algorithm::JetAlgorithm.Algorithm)

Check if the algorithm is a e+e- reconstruction algorithm.

# Returns
`true` if the algorithm is a e+e- reconstruction algorithm, `false` otherwise.
"""
function is_ee(algorithm::JetAlgorithm.Algorithm)
    return algorithm in (JetAlgorithm.EEKt, JetAlgorithm.Durham, JetAlgorithm.Valencia)
end

"""
    enum RecombinationScheme

An EnumX scoped enumeration representing different recombination schemes that
are supported directly in the package.

These schemes map to both a `recombine` and a `preprocess` function, which are
used in the main reconstruction algorithm.
"""
@enumx T=Recombine RecombinationScheme ESchemeRaw EScheme PtScheme Pt2Scheme
const AllRecombinationSchemes = [String(Symbol(x))
                                 for x in instances(RecombinationScheme.Recombine)]

# Note it's a bit fragile to have the dictionary and the enum built
# separately, but it is manageable. There is a test in the CI that
# checks that all the enums are defined in the dictionary.
const RecombinationMethods = Dict(RecombinationScheme.EScheme => (recombine = addjets_escheme,
                                                                  preprocess = preprocess_escheme),
                                  RecombinationScheme.ESchemeRaw => (recombine = addjets_escheme,
                                                                     preprocess = nothing),
                                  RecombinationScheme.PtScheme => (recombine = addjets_ptscheme,
                                                                   preprocess = preprocess_ptscheme),
                                  RecombinationScheme.Pt2Scheme => (recombine = addjets_pt2scheme,
                                                                    preprocess = preprocess_pt2scheme))
