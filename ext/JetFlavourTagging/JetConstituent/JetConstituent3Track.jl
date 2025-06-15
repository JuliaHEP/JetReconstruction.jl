using EDM4hep
using JetReconstruction
using StructArrays: StructVector


### Track Related Functions (10)

## Track Parameters
# get_d0 - Impact parameter in transverse plane
# get_z0 - Impact parameter along z-axis
# get_phi0 - Initial azimuthal angle
# get_omega - Track curvature
# get_tanLambda - Tangent of dip angle

## Track Parameter Transformations (XPtoPar)
# XPtoPar_dxy - Transformed transverse impact parameter
# XPtoPar_dz - Transformed longitudinal impact parameter
# XPtoPar_phi - Transformed azimuthal angle
# XPtoPar_C - Track curvature parameter
# XPtoPar_ct - c×tau parameter