# Jet Substructure

## Structures

### `MassDropTagger`

The `MassDropTagger` structure is used for tagging jets based on mass drop conditions, which helps in identifying subjets within a jet that undergo a significant drop in mass.

**Fields**:

- `mu::Float64`: Maximum allowed mass ratio for a jet to pass the tagging condition.
- `y::Float64`: Minimum kT distance threshold for parent jet separation.

---

### `SoftDropTagger`

The `SoftDropTagger` instance is used to apply soft-drop grooming to jets, removing soft, wide-angle radiation. This approach is commonly used in jet grooming to reduce contamination from soft particles.

**Fields**:

- `zcut::Float64`: Minimum allowed energy fraction for subjets.
- `b::Float64`: Angular exponent controlling soft radiation suppression.
- `cluster_rad::Float64`: New radius used to recluster components of the jet. Defaults to `1.0` if no value is specified.

---

### `JetFilter`

The `JetFilter` structure is used to filter jets based on a specific radius and the number of hardest subjets. This technique reduces contamination from peripheral soft particles.

**Fields**:

- `filter_radius::Float64`: Radius parameter used to recluster subjets.
- `num_hardest_jets::Int`: Number of hardest subjets retained in the filtered result.

---

### `JetTrim`

`JetTrim` instance is used to trim jets by removing soft, large-angle components from the jet. This is useful in cleaning up jets to remove softer particles at wide angles.

**Fields**:

- `trim_radius::Float64`: Radius used for reclustering in trimming.
- `trim_fraction::Float64`: Minimum momentum fraction for retained subjets.
- `recluster_method::JetAlgorithm.Algorithm`: Method identifier for reclustering.

---

## Functions

### `mass_drop`

```julia
mass_drop(jet::PseudoJet, clusterseq::ClusterSequence, tag::MassDropTagger) -> PseudoJet
```

The `mass_drop` function identifies subjets in a jet that pass the mass drop tagging condition. It iterates through the clustering history of the jet, stopping at the first jet that satisfies the mass and distance thresholds.

**Arguments** :

* `jet`: `PseudoJet` instance representing the jet to be tagged.
* `clusterseq`: `ClusterSequence` with jet clustering history.
* `tag`: `MassDropTagger` instance providing mass drop parameters.

**Returns** :

`PseudoJet`: The jet (or subjet) that satisfies the mass drop condition, or a zero-momentum `PseudoJet` if no tagging occurs.

---

### `soft_drop`

```julia
soft_drop(jet::PseudoJet, clusterseq::ClusterSequence, tag::SoftDropTagger) -> PseudoJet
```

The `soft_drop` function applies soft-drop grooming to remove soft, wide-angle radiation from jets. It reclusters the jet with a specified radius and clustering method, iteratively checking the soft-drop condition on subjets.

**Arguments** :

* `jet`: `PseudoJet` instance to groom.
* `clusterseq`: `ClusterSequence` containing jet history.
* `tag`: `SoftDropTagger` instance with soft-drop parameters.

**Returns** :

`PseudoJet`: Groomed jet or zero-momentum `PseudoJet` if grooming fails.

---

### `jet_filtering`

```julia-repl
jet_filtering(jet::PseudoJet, clusterseq::ClusterSequence, filter::JetFilter) -> PseudoJet
```

The `jet_filtering` function filters a jet to retain only the hardest subjets based on a specified radius and number. This helps in refining the jet structure by reducing soft particle contamination.

**Arguments** :

* `jet`: `PseudoJet` instance representing the jet to filter.
* `clusterseq`: `ClusterSequence` containing jet history.
* `filter`: `JetFilter` instance specifying radius and number of subjets.

**Returns** :

`PseudoJet`: Filtered jet composed of the hardest subjets.

---

### `jet_trimming`

```julia
jet_trimming(jet::PseudoJet, clusterseq::ClusterSequence, trim::JetTrim) -> PseudoJet
```

The `jet_trimming` function trims a jet by removing subjets with transverse momentum below a specified fraction of the main jet's momentum. This method cleans up jets by removing soft particles.

**Arguments** :

* `jet`: `PseudoJet` instance representing the jet to trim.
* `clusterseq`: `ClusterSequence` containing jet history.
* `trim`: `JetTrim` instance specifying trimming parameters.

**Returns** :

`PseudoJet`: Trimmed jet composed of retained subjets.
