"""
    jet_reconstruct(particles::AbstractVector; algorithm::JetAlgorithm.Algorithm,
                         p::Union{Real, Nothing} = nothing, R = 1.0,
                         recombine = addjets, preprocess = nothing,
                         strategy::RecoStrategy.Strategy = RecoStrategy.Best)

Reconstructs jets from a collection of particles using a specified algorithm and
strategy.

# Arguments
- `particles::AbstractVector`: A collection of particles used for jet
  reconstruction. 
- `algorithm::JetAlgorithm.Algorithm`: The algorithm to use for jet reconstruction.
- `p::Union{Real, Nothing} = nothing`: The power value used for the distance
  measure for generalised k_T algorithms (GenKt, EEKt). Other algorithms will
  ignore this value.
- `R = 1.0`: The jet radius parameter.
- `recombine = addjets`: The recombination scheme used for combining particles.
- `preprocess = nothing`: The function to preprocess the particles before
  reconstruction (e.g., for massless schemes). `nothing` means the particles are
  not preprocessed.
- `strategy::RecoStrategy.Strategy = RecoStrategy.Best`: The jet reconstruction
   strategy to use. `RecoStrategy.Best` makes a dynamic decision based on the
   number of starting particles.

Note that `p` must be specified for `GenKt` and `EEKt` algorithms,
other algorithms will ignore its value.
When an algorithm has no `R` dependence the `R` parameter is ignored.

# Returns
A cluster sequence object containing the reconstructed jets and the merging
history.

# Details

## `particles` argument
Any type that supplies the methods `pt2()`, `phi()`, `rapidity()`, `px()`,
`py()`, `pz()`, `energy()` (in the `JetReconstruction` namespace) can be used.
This includes `LorentzVector`, `LorentzVectorCyl`, `PseudoJet` and `EEJet`,
for which these methods are already predefined in the `JetReconstruction` namespace.

**Note** when using `PseudoJet` or `EEJet`, the history indices (`_cluster_hist_index`)
must be set correctly. For initial jets, this means assigning to each jet its index
in the vector.
    
## `recombine` argument
The `recombine` argument is the function used to merge pairs of particles. The
default is `addjets`, which uses 4-momenta addition (a.k.a. the E-scheme).

## `preprocess` argument
The `preprocess` argument is a function that will be called for all original
input particles and which returns a new particle, usually matching a
non-standard recombination scheme, e.g., massless particles for ``p_T`` or
``p_T^2`` recombination. `nothing` means no preprocessing is done.

# Example
```julia
jet_reconstruct(particles; algorithm = JetAlgorithm.AntiKt, R = 0.4)
jet_reconstruct(particles; algorithm = JetAlgorithm.Kt, R = 1.0)
jet_reconstruct(particles; algorithm = JetAlgorithm.Durham)
jet_reconstruct(particles; algorithm = JetAlgorithm.GenKt, p = 0.5, R = 1.0)
jet_reconstruct(particles; algorithm = JetAlgorithm.AntiKt, R = 1.0, preprocess = preprocess_ptscheme, 
                recombine = addjets_ptscheme)
```
"""
function jet_reconstruct(particles::AbstractVector; algorithm::JetAlgorithm.Algorithm,
                         p::Union{Real, Nothing} = nothing, R = 1.0,
                         recombine = addjets, preprocess = nothing,
                         strategy::RecoStrategy.Strategy = RecoStrategy.Best)
    if is_pp(algorithm)
        # We assume a pp reconstruction
        if strategy == RecoStrategy.Best
            # The breakpoint of ~80 is determined empirically on e+e- -> H and 0.5 TeV pp -> 5 GeV jets
            alg = length(particles) > 80 ? tiled_jet_reconstruct : plain_jet_reconstruct
        elseif strategy == RecoStrategy.N2Plain
            alg = plain_jet_reconstruct
        elseif strategy == RecoStrategy.N2Tiled
            alg = tiled_jet_reconstruct
        else
            throw(ErrorException("Invalid strategy: $(strategy)"))
        end
    elseif is_ee(algorithm)
        alg = ee_genkt_algorithm
    else
        throw(ErrorException("Invalid algorithm neither pp nor ee: $(algorithm)"))
    end

    # Now call the chosen algorithm, passing through the other parameters
    alg(particles; p = p, algorithm = algorithm, R = R, recombine = recombine,
        preprocess = preprocess)
end
