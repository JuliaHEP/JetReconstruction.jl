using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### Particle Identification (7)

# get_PIDs - Get particle IDs from MC matching
# get_PIDs_cluster - For clustered jets
# get_isMu - Check if constituent is muon
# get_isEl - Check if constituent is electron
# get_isChargedHad - Check if constituent is charged hadron
# get_isGamma - Check if constituent is photon
# get_isNeutralHad - Check if constituent is neutral hadron