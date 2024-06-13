"""
    jet_reconstruct(particles; p = -1, R = 1.0, recombine = +, strategy = RecoStrategy.Best)

Reconstructs jets from a collection of particles using a specified algorithm and
strategy

# Arguments
- `particles`: A collection of particles used for jet reconstruction. 
- `p=-1`: The power value used for the distance measure for generalised k_T,
  which maps to a particular reconstruction algorithm (-1 = AntiKt, 0 =
  Cambridge/Aachen, 1 = Kt).
- `R=1.0`: The jet radius parameter.
- `recombine=+`: The recombination scheme used for combining particles.
- `strategy=RecoStrategy.Best`: The jet reconstruction strategy to use.
  `RecoStrategy.Best` makes a dynamic decision based on the number of starting
  particles.

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

# Example
```julia
jet_reconstruct(particles; p = -1, R = 0.4)
```
"""
function jet_reconstruct(particles; p = -1, R = 1.0, recombine = +,
                         strategy = RecoStrategy.Best)
    # Either map to the fixed algorithm corresponding to the strategy
    # or to an optimal choice based on the density of initial particles

    if strategy == RecoStrategy.Best
        # The breakpoint of ~80 is determined empirically on e+e- -> H and 0.5TeV pp -> 5GeV jets
        algorithm = length(particles) > 80 ? tiled_jet_reconstruct : plain_jet_reconstruct
    elseif strategy == RecoStrategy.N2Plain
        algorithm = plain_jet_reconstruct
    elseif strategy == RecoStrategy.N2Tiled
        algorithm = tiled_jet_reconstruct
    else
        throw(ErrorException("Invalid strategy: $(strategy)"))
    end

    # Now call the chosen algorithm, passing through the other parameters
    algorithm(particles; p = p, R = R, recombine = recombine)
end
