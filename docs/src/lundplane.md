# Lund Jet Plane

The Lund jet plane ([arXiv:1807.04758](https://arxiv.org/abs/1807.04758)) is a simple and intuitive way to visualize how a jet breaks up into smaller components. Each point in the plane represents a branching in the jet, helping understand how energy and angles are distributed within.

---

### Decluster

```julia
decluster(jet::T, clusterseq::ClusterSequence{T}) where {T <: FourMomentum} -> Tuple{PseudoJet}
```

Recursively declusters a jet into two parent subjets using a provided clustering sequence. The subjets are ordered by decreasing transverse momentum squared (`pt²`).

```julia
j1, j2 = decluster(jet, clusterseq)
```

- If no parents are found, the function returns `nothing` for one or both subjets.

Useful for building recursive jet trees and jet grooming algorithms.

---

### Generate Lund Emissions

```julia
generate_lund_emissions(jet::PseudoJet, cs::ClusterSequence{PseudoJet}) -> Vector{NamedTuple}
```

Constructs the Lund plane emissions of a jet. The result is a list of declustering steps, each represented by a named tuple of physics observables.

```julia
lund_points = generate_lund_emissions(jet, cluster_seq)
```

Each named tuple includes:

- `h_pt`: transverse momentum of the harder branch
- `s_pt`: transverse momentum of the softer branch
- `z`: momentum fraction
- `delta`: angular distance ΔR
- `kt`: transverse momentum of the softer branch relative to the harder
- `psi`: azimuthal angle between the two splittings
- `kappa`: z × ΔR

These values are useful for visualization and physics analysis of jet splittings.

---

### Visualisation of Lund Jet Plane

`examples/lundplane/lund-plane-visualisation.jl` provides an example for generating and plotting the average lund image.

```julia
generate_average_lund_image(njets::Int, delta_array::Vector{Vector{Real}}, kt_array::Vector{Vector{Real}}; xrange::Tuple{Real}, yrange::Tuple{Real}, bins::Int) -> (xgrid, ygrid, avg_image)
```

Computes the average Lund image over a set of jets, by binning (log(1/ΔR), log(kt)) values from each jet into a 2D histogram and averaging the histograms per jet.

```julia
x, y, avg_image = generate_average_lund_image(njets, delta_array, kt_array; bins=25)
```

- `delta_array` and `kt_array` should be arrays of arrays: one per jet.
- The output is a tuple of x and y bin edges, and a 2D image array.
- This is useful for creating heatmaps to compare jet substructure distributions.

---

For further examples and advanced use cases, see the documentation in the `examples/lundplane/` directory
