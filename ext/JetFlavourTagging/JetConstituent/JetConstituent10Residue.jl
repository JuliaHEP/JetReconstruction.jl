using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### Four-Vector Operations and Residuals (11)

# compute_tlv_jets - Convert jets to TLorentzVector
# sum_tlv_constituents - Sum constituent four-vectors
# InvariantMass - Calculate invariant mass of two objects
# all_invariant_masses - All pairwise invariant masses
# compute_residue_energy - Energy difference (constituents - jet)
# compute_residue_pt - pT difference
# compute_residue_phi - φ difference
# compute_residue_theta - θ difference
# compute_residue_px - px difference
# compute_residue_py - py difference
# compute_residue_pz - pz difference