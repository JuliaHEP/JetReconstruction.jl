## Jet class constructors from EMD4hep objects

module EDM4hepJets

using JetReconstruction
using EDM4hep

"""
    JetReconstruction.px(recoparticle::ReconstructedParticle)

Return the x component of the momentum of a ReconstructedParticle.
"""
JetReconstruction.px(recoparticle::ReconstructedParticle) = recoparticle.momentum.x

"""
    JetReconstruction.py(recoparticle::ReconstructedParticle)

Return the y component of the momentum of a ReconstructedParticle.
"""
JetReconstruction.py(recoparticle::ReconstructedParticle) = recoparticle.momentum.y

"""
    JetReconstruction.pz(recoparticle::ReconstructedParticle)

Return the z component of the momentum of a ReconstructedParticle.
"""
JetReconstruction.pz(recoparticle::ReconstructedParticle) = recoparticle.momentum.z

"""
    JetReconstruction.energy(recoparticle::ReconstructedParticle)

Return the energy component of a ReconstructedParticle's four vector.
"""
JetReconstruction.energy(recoparticle::ReconstructedParticle) = recoparticle.energy

"""
    JetReconstruction.EEJet(recoparticle::ReconstructedParticle)

Construct an EEJet from a ReconstructedParticle.
"""
function JetReconstruction.EEJet(recoparticle::ReconstructedParticle)
    EEJet(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
          JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle))
end

"""
    JetReconstruction.PseudoJet(recoparticle::ReconstructedParticle)

Construct an PseudoJet from a ReconstructedParticle.
"""
function JetReconstruction.PseudoJet(recoparticle::ReconstructedParticle)
    PseudoJet(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
              JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle))
end

end
