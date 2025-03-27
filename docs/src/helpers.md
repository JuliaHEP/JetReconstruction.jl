# Jet Helper Functions

These functions are provided as a convenient way to work with the results of jet
reconstruction.

## Jet Pair Helpers

- [`pt_fraction`](@ref)

Returns the transverse momentum fraction in the softer of the two jets.

- [`kt_scale`](@ref)

Returns the transverse momentum scale between two jets, the product of the
smaller ``p_T`` and the angular separation in the the η-ϕ plane.

## Conversion Functions

- [`lorentzvector`](@ref)

Return a cartesian `LorentzVector` from a jet.

- [`lorentzvector_cyl`](@ref)

Return a cylindrical `LorentzVectorCyl` from a jet.
