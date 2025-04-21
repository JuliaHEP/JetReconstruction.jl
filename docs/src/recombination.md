# Jet Recombination

When two jets are merged different strategies can be adopted to produce the merged jet.

The function to do this can be set by the user and passed as a parameter to the
reconstruction algorithms.

## Default - Four Vector Addition

The default for jet merging is simply four momentum addition, that is:

``
(\mathbf{p}_m, E_m) = (\mathbf{p_1} + \mathbf{p_2}, E_1 + E_2)
``

This is defined as the [`addjets`](@ref) function in the package, which also
serves as an example of how these functions are written.

### Different Recombination Schemes

Two additional recombination schemes are supported, the ``p_T`` and ``p_T^2``
schemes. In these schemes the recombined jet is created to be *massless*, i.e.,
the mass is set to the 3-momentum. The transverse momentum is the sum of the two
parent jets and the rapidity (``y``) and phi (``\phi``) values are weighted
averages, by ``p_T`` or ``p_T^2``, of the parent jets.

- [`addjets_ptscheme`](@ref)
- [`addjets_pt2scheme`](@ref)

## User Defined Recombination

If a different merging scheme is desired then a method must be defined
that implements the following interface:

```julia
user_recombine(jet1::T, jet2::T, cluster_hist_index::Int) where {T <: FourMomentum} -> T
```

i.e., three arguments are needed, the two parent jets and the
`cluster_hist_index`, which is needed to identify the jet in the reconstruction
sequence.

It is recommended to use the constructor signature for the output jet of:

```julia
T(px, py, pz, E, cluster_hist_index)
```

Where `px`, `py`, `pz` and `E` have been calculated from the inputs `jet1` and
`jet2` as desired.

However, if working in ``(p_T, y, \phi, m)`` space, use the alternative constructor
with named parameters:

```julia
T(;pt=pt, rap=rap, phi=phi, m=m, cluster_hist_index=cluster_hist_index)
```

(Note that there is a default of `m=0.0`, which is used for massless
recombination.)

The user function should not modify the `cluster_hist_index`, but must pass in
to the new jet's constructor to ensure that the resulting reconstruction
[`ClusterSequence`](@ref) is valid. The recombination functions defined in the
package serve as examples.

### Using an Custom Recombination Method

To use a non-default recombination method, simply pass the method to the
[`jet_reconstruct`](@ref) entry point as the `recombine` parameter.
