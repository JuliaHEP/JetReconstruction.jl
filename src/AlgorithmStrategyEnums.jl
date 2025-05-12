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
"""
@enumx T=Algorithm JetAlgorithm AntiKt CA Kt GenKt EEKt Durham
const AllJetRecoAlgorithms = [String(Symbol(x)) for x in instances(JetAlgorithm.Algorithm)]

"""
    const varpower_algorithms

A constant array that contains the jet algorithms for which power is variable.

"""
const varpower_algorithms = [JetAlgorithm.GenKt, JetAlgorithm.EEKt]

"""
    power2algorithm

A dictionary that maps power values to corresponding jet algorithm used for pp
jet reconstruction, where these are fixed.
"""
const power2algorithm = Dict(-1 => JetAlgorithm.AntiKt,
                             0 => JetAlgorithm.CA,
                             1 => JetAlgorithm.Kt)

"""
    algorithm2power

A dictionary that maps pp algorithm names to their corresponding power values.

The dictionary is created by iterating over the `power2algorithm` dictionary and swapping the keys and values.
"""
const algorithm2power = Dict((i.second, i.first) for i in power2algorithm)

"""
    get_algorithm_power_consistency(; p::Union{Real, Nothing}, algorithm::Union{JetAlgorithm.Algorithm, Nothing})

Get the algorithm and power consistency correct

This function checks the consistency between the algorithm and power parameters.
If the algorithm is specified, it checks if the power parameter is consistent
with the algorithm's known power. If the power parameter is not specified, it
sets the power parameter based on the algorithm. If neither the algorithm nor
the power parameter is specified, it throws an `ArgumentError`.

# Arguments
- `p::Union{Real, Nothing}`: The power value.
- `algorithm::Union{JetAlgorithm.Algorithm, Nothing}`: The algorithm.

# Returns
A named tuple of the consistent power and algorithm values.

# Throws
- `ArgumentError`: If the algorithm and power are inconsistent or if neither the
  algorithm nor the power is specified.

"""
function get_algorithm_power_consistency(; p::Union{Real, Nothing},
                                         algorithm::Union{JetAlgorithm.Algorithm, Nothing})
    # For the case where an algorithm has a variable power value
    # we need to check that the power value and algorithm are both specified
    if algorithm in varpower_algorithms
        if isnothing(p)
            throw(ArgumentError("Power must be specified for algorithm $algorithm"))
        end
        # And then we just return what these are...
        return (p = p, algorithm = algorithm)
    end

    # Some algorithms have fixed power values
    if algorithm == JetAlgorithm.Durham
        return (p = 1, algorithm = algorithm)
    end

    # Otherwise we check the consistency between the algorithm and power
    if !isnothing(algorithm)
        power_from_alg = algorithm2power[algorithm]
        if !isnothing(p) && p != power_from_alg
            throw(ArgumentError("Algorithm and power are inconsistent"))
        end
        return (p = power_from_alg, algorithm = algorithm)
    else
        if isnothing(p)
            throw(ArgumentError("Either algorithm or power must be specified"))
        end
        # Set the algorithm from the power
        algorithm_from_power = power2algorithm[p]
        return (p = p, algorithm = algorithm_from_power)
    end
end

"""Allow a check for algorithm and power consistency"""
function check_algorithm_power_consistency(; p::Union{Real, Nothing},
                                           algorithm::Union{JetAlgorithm.Algorithm,
                                                            Nothing})
    get_algorithm_power_consistency(p = p, algorithm = algorithm)
    return true
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
    return algorithm in [JetAlgorithm.EEKt, JetAlgorithm.Durham]
end

"""
    enum RecombinationScheme

An EnumX scoped enumeration representing different recombination schemes that
are supported directly in the package.

These schemes map to both a `recombine` and a `preprocess` function, which are
used in the main reconstruction algorithm.
"""
@enumx T=Recombine RecombinationScheme EScheme PtScheme Pt2Scheme
const AllRecombinationSchemes = [String(Symbol(x))
                                 for x in instances(RecombinationScheme.Recombine)]

const RecombinationMethods = Dict(RecombinationScheme.EScheme => (recombine = addjets_escheme,
                                                                  preprocess = nothing),
                                  RecombinationScheme.PtScheme => (recombine = addjets_ptscheme,
                                                                   preprocess = preprocess_ptscheme),
                                  RecombinationScheme.Pt2Scheme => (recombine = addjets_pt2scheme,
                                                                    preprocess = preprocess_pt2scheme))
