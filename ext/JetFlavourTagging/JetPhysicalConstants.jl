"""
Physical constants and special values used in the JetFlavourTagging extension.

This module contains all physical constants, particle properties, tolerance values,
and special markers used throughout the jet flavour tagging algorithms.
"""
module JetPhysicalConstants

using PhysicalConstants.CODATA2018: c_0, m_e, ħ, e, k_B, N_A, R

# Note: PhysicalConstants.jl provides values with units and uncertainties
# To extract just the numerical value, use .val
# Example: c_0.val gives the speed of light value without units

# Physical Constants from PhysicalConstants.jl
# Extract the numerical value and convert to Float32
const C_LIGHT = Float32(c_0.val)  # Speed of light in m/s from CODATA2018
const C_LIGHT_INV = 1.0f0 / C_LIGHT  # Inverse speed of light

# Particle Masses (in GeV/c²)
# Electron mass from CODATA2018 (convert from kg to GeV/c²)
# 1 GeV/c² = 1.78266192e-27 kg
const ELECTRON_MASS = Float32(m_e.val / 1.78266192e-27)  # Electron mass in GeV/c²
# Note: Muon and pion masses are not in CODATA2018, using PDG values
const MUON_MASS = 0.105658f0  # Muon mass in GeV/c² (PDG value)
const PION_MASS = 0.13957039f0  # Charged pion mass (π±) in GeV/c² (PDG value)

# Mass Comparison Tolerances
const ELECTRON_TOLERANCE = 1.0f-5  # Tolerance for electron mass comparison
const MUON_TOLERANCE = 1.0f-3  # Tolerance for muon mass comparison
const PION_TOLERANCE = 1.0f-3  # Tolerance for pion mass comparison

# PDG Particle ID Codes
const PDG_PHOTON = 22  # Photon (γ)
const PDG_K_LONG = 130  # K⁰_L (K-long neutral hadron)

# Special/Undefined Values
const UNDEF_VAL = -9.0f0  # Sentinel value for missing/invalid data
const INVALID_TOF_MASS = 9.0f0  # Invalid mass from time-of-flight
const INVALID_MASS = -1.0f0  # Invalid mass calculation marker

# Unit Conversion Factors
const MM_TO_M = 0.001f0  # Millimeter to meter conversion
const NS_TO_S = 1.0f-9  # Nanosecond to second conversion
const PS_TO_S = 1.0f-12  # Picosecond to second conversion
const FS_TO_S = 1.0f-15  # Femtosecond to second conversion

# Export all constants
export C_LIGHT, C_LIGHT_INV
export ELECTRON_MASS, MUON_MASS, PION_MASS
export ELECTRON_TOLERANCE, MUON_TOLERANCE, PION_TOLERANCE
export PDG_PHOTON, PDG_K_LONG
export UNDEF_VAL, INVALID_TOF_MASS, INVALID_MASS
export MM_TO_M, NS_TO_S, PS_TO_S, FS_TO_S

end # module JetPhysicalConstants