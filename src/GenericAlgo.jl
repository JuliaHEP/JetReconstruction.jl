"""
    jet_reconstruct(particles; p::Union{Real, Nothing} = nothing,
                         algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                         R = 1.0, recombine = +,
                         strategy::RecoStrategy.Strategy = RecoStrategy.Best)

Reconstructs jets from a collection of particles using a specified algorithm and
strategy.

# Arguments
- `particles`: A collection of particles used for jet reconstruction. 
- `p::Union{Real, Nothing} = nothing`: The power value used for the distance measure
  for generalised k_T, which maps to a particular reconstruction algorithm (-1 =
  AntiKt, 0 = Cambridge/Aachen, 1 = Kt).
- `algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing`: The algorithm
  to use for jet reconstruction.
- `R = 1.0`: The jet radius parameter.
- `recombine = +`: The recombination scheme used for combining particles.
- `strategy::RecoStrategy.Strategy = RecoStrategy.Best`: The jet reconstruction
   strategy to use. `RecoStrategy.Best` makes a dynamic decision based on the
   number of starting particles.

Note that one of `p` or `algorithm` must be specified, with `algorithm` preferred.

# Returns
A cluster sequence object containing the reconstructed jets and the merging
history.

# Details

## `particles` argument
Any type that supplies the methods `pt2()`, `phi()`, `rapidity()`, `px()`,
`py()`, `pz()`, `energy()` (in the `JetReconstruction` namespace) can be used.
This includes `LorentzVector`, `LorentzVectorCyl`, and `PseudoJet`, for which
these methods are already predefined in the `JetReconstruction` namespace.

## `recombine` argument
The `recombine` argument is the function used to merge pairs of particles. The
default is simply `+(jet1,jet2)`, i.e. 4-momenta addition or the *E*-scheme.

## Consistency of `p`, `algorithm` and `R` arguments
If an algorithm is explicitly specified the `p` value should be consistent with
it or `nothing`. If the algorithm is one where `p` can vary, then it has to be
given, along with the algorithm.

If the `p` parameter is passed and `algorithm=nothing`, then pp-type
reconstruction is implied (i.e., AntiKt, CA, Kt or GenKt will be used,
depending on the value of `p`).

When an algorithm has no `R` dependence the `R` parameter is ignored.

# Example
```julia
jet_reconstruct(particles; p = -1, R = 0.4)
jet_reconstruct(particles; algorithm = JetAlgorithm.Kt, R = 1.0)
jet_reconstruct(particles; algorithm = JetAlgorithm.Durham)
jet_reconstruct(particles; algorithm = JetAlgorithm.GenKt, p = 0.5, R = 1.0)
```
"""
function jet_reconstruct(particles; p::Union{Real, Nothing} = nothing,
                         algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                         R = 1.0, recombine = +,
                         strategy::RecoStrategy.Strategy = RecoStrategy.Best)

    # Either map to the fixed algorithm corresponding to the strategy
    # or to an optimal choice based on the density of initial particles
    if (algorithm === nothing) && (p === nothing)
        throw(ArgumentError("Please specify either an algorithm or a power value (power only for pp-type reconstruction)"))
    end

    if (algorithm === nothing) || is_pp(algorithm)

        # We assume a pp reconstruction
        if strategy == RecoStrategy.Best
            # The breakpoint of ~80 is determined empirically on e+e- -> H and 0.5 TeV pp -> 5 GeV jets
            alg = length(particles) > 800 ? tiled_jet_reconstruct : plain_jet_reconstruct
        elseif strategy == RecoStrategy.N2Plain
            alg = plain_jet_reconstruct
        elseif strategy == RecoStrategy.N2Tiled
            alg = tiled_jet_reconstruct
        else
            throw(ErrorException("Inv alid strategy: $(strategy)"))
        end
    elseif is_ee(algorithm)
        alg = ee_genkt_algorithm
    end

    # Now call the chosen algorithm, passing through the other parameters
    alg(particles; p = p, algorithm = algorithm, R = R, recombine = recombine)

end
