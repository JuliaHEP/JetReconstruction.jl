# Jet Recombination

When two jets are merged different strategies can be adopted to produce the merged jet.

There are two functions to support this, which can be set by the user and passed
as a parameters to the reconstruction algorithms. One function gives the
necessary *preprocessing* for input particles, e.g., setting particles to be
massless. The other controls the actual *recombination* of two particles into a
merged jet.

These functions are passed as the `preprocess` and `recombine` parameters to the
reconstruction interfaces.

## Default - Four Vector Addition

The default for jet merging is simply four momentum addition, that is:

``
(\mathbf{p}_m, E_m) = (\mathbf{p_1} + \mathbf{p_2}, E_1 + E_2)
``

This is defined as the [`addjets`](@ref) function in the package, which also
serves as an example of how the recombination functions are written.

In this case, the preprocessing [`preprocess_escheme`](@ref) is limited to assigning history index to the jet. Should the jets have the correct history index already, no preprocessing is needed and `preprocess = nothing` can be used.

### Different Recombination Schemes

Two additional recombination schemes are directly supported, the ``p_T`` and
``p_T^2`` schemes. In these schemes the recombined jet is created to be
*massless*, i.e., the mass is set to the 3-momentum. The transverse momentum is
the sum of the two parent jets and the rapidity (``y``) and phi (``\phi``)
values are weighted averages, by ``p_T`` or ``p_T^2``, of the parent jets.

- `recombine =` [`addjets_ptscheme`](@ref)
- `recombine =` [`addjets_pt2scheme`](@ref)

In this case the input particles must be rescaled to be massless, setting the
energy equal to the (three) momentum sum.

- `preprocess =` [`preprocess_ptscheme`](@ref)
- `preprocess =` [`preprocess_pt2scheme`](@ref)

(In fact `preprocess_pt2scheme` is just an alias for `preprocess_ptscheme` as
the rescaling is identical.)

### Named Recombination Schemes

To simplify the usage of different recombination schemes supported directly,
there is a defined enum (scoped, using `EnumX`) for each one:
`RecombinationScheme.SCHEME`.

This enum is then used with the `RecombinationMethods` dictionary to
obtain a named tuple in which `recombine` and `preprocess` are set, which can
then be splatted into the [`jet_reconstruct`](@ref) interface:

```julia
myscheme = RecombinationMethods[RecombinationScheme.PtScheme]
jet_reconstruct(event; R = distance, p = p, algorithm = algorithm,
                                 strategy = strategy, myscheme...)
```

The supported values in the enum are:

| Scheme | Implements |
|---|---|
| `EScheme` | Default 4-momentum addition |
| `ESchemeRaw` | Four-momentum addition without preprocessing history index |
| `PtScheme` | Massless weighted average of momentum |
| `Pt2Scheme` | Massless weighted average of momentum squared |

All schemes except for `ESchemeRaw` include the necessary preprocessing to set the
history index of the input particles.
(Should other schemes prove to be particularly desired they can be implemented
on request.)

## User Defined Recombination

### Preprocessing

The user must supply, if needed, a preprocessing function, which accepts an
input particle and output type, and returns the rescaled particle. This function must accept a named argument `cluster_hist_index` to pass to the constructor of the resulting
particle. `EEJet` or `PseudoJet` should be valid output types depending whether ``e^+e^-`` or ``pp`` reconstruction is being performed.

```julia
user_preprocess(jet::T, ::Type{OutputT}; cluster_hist_index) -> OutputT
```

An example of a preprocessing function is [`preprocess_ptscheme`](@ref).

### Recombination

If a different merging scheme is desired then a method must be defined
that implements the following interface:

```julia
user_recombine(jet1::T, jet2::T; cluster_hist_index::Int) where {T <: FourMomentum} -> T
```

i.e., three arguments are needed, the two parent jets and the named argument
`cluster_hist_index`, which is needed to identify the jet in the reconstruction
sequence.

It is recommended to use the constructor signature for the output jet of:

```julia
T(px, py, pz, E; cluster_hist_index = cluster_hist_index)
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
package serve as examples: [`addjets_ptscheme`](@ref).

### Using an Custom Recombination Method

To use a non-default recombination method, simply pass the recombination method
to the [`jet_reconstruct`](@ref) entry point as the `recombine` parameter and
the preprocessing method as `preprocess`.

A very convenient way to do this is to bind these functions into a named tuple
and splat the tuple into the arguments for the reconstruction.
