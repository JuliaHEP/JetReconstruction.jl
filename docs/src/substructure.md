# Jet Substructure

Jet substructure techniques provide powerful tools for analyzing and refining the properties of jets. Below are some of the key jet substructure functions, that are available.

### Mass Drop Tagging

```julia
mass_drop(jet::PseudoJet, clusterseq::ClusterSequence; mu::Float64, y::Float64) -> PseudoJet
```

The [mass_drop](@ref) function identifies subjets in a jet that pass the mass drop tagging condition. To use the `mass_drop` function:

- Apply the `mass_drop` function to the jet and its clustering sequence.

```julia
tagger = (mu=0.67, y=0.09)
tagged_jet = mass_drop(jet, clusterseq; tagger...)
```
or
```julia
tagged_jet = mass_drop(jet, clusterseq; mu = 0.67, y = 0.09)
```

- If the jet is tagged successfully, the function returns the identified subjet. Else it returns `PseudoJet(0.0, 0.0, 0.0, 0.0)`.

---

### Soft Drop Tagging

```julia
soft_drop(jet::PseudoJet, clusterseq::ClusterSequence, tag::SoftDropTagger) -> PseudoJet
```

The [soft_drop](@ref) function applies soft-drop grooming to remove soft, wide-angle radiation from jets. It reclusters the jet with a specified radius and clustering method, iteratively checking the soft-drop condition on subjets. To use the `soft_drop` function:

- Apply the `soft_drop` function to the jet and its clustering sequence.

```julia
tagger = (zcut = 0.1, beta = 2.0)
tagged_jet = soft_drop(jet, clusterseq; tagger...)
```
or
```julia
tagged_jet = soft_drop(jet, clusterseq; zcut = 0.1, beta = 2.0)
```

By default, the reclustering radius is set to `1.0` which can be modified by the user as:
```julia
tagger = (zcut = 0.1, beta = 2.0, radius = 0.4)
tagged_jet = soft_drop(jet, clusterseq; tagger...)
```
or
```julia
tagged_jet = soft_drop(jet, clusterseq; zcut = 0.1, beta = 2.0, radius = 0.4)
```

- If the jet is tagged successfully, the function returns the identified subjet. Else it returns `nothing`.

---

### Jet Filtering

```julia
jet_filtering(jet::PseudoJet, clusterseq::ClusterSequence, filter::JetFilter) -> PseudoJet
```

The [jet_filtering](@ref) function filters a jet to retain only the hardest subjets based on a specified radius and number. To use the function:

- Apply the `jet_filtering` function to refine the jet.

```julia
filter = (radius=0.3, hardest_jets=3)
filtered_jet = jet_filtering(jet, clusterseq; filter...)
```
or
```julia
filtered_jet = jet_filtering(jet, clusterseq; radius = 0.3, hardest_jets = 3)
```

- The function returns the filtered jet that retains only the most significant subjets, reducing noise.

---

### Jet Trimming

```julia
jet_trimming(jet::PseudoJet, clusterseq::ClusterSequence, trim::JetTrim) -> PseudoJet
```

The [jet_trimming](@ref) function trims a jet by removing subjets with transverse momentum below a specified fraction of the main jet's momentum. This method cleans up jets by removing soft particles. To use this function:

- Apply the `jet_trimming` function to clean the jet.

```julia
trim = (radius=0.3, fraction=0.3, recluster_method=JetAlgorithm.CA)
trimmed_jet = jet_trimming(jet, clusterseq; trim...)
```
or
```julia
trimmed_jet = jet_trimming(jet, clusterseq; radius=0.3, fraction=0.3, recluster_method=JetAlgorithm.CA)
```

- The function returns the trimmed jet.

It is to be noted that the `jet_trimming` function reclusters the constituents of the jet using either `C/A` or `kT ` algorithm, which needs to be specified.

---

For a more detailed guide to use these substructure modules, one can refer to the provided examples in the `examples/substructure` directory.
