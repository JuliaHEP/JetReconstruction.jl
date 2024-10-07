## Jet class constructors from EMD4hep objects

module EDM4hepJets

using JetReconstruction
using EDM4hep

JetReconstruction.px(recoparticle::ReconstructedParticle) = recoparticle.momentum.x
JetReconstruction.py(recoparticle::ReconstructedParticle) = recoparticle.momentum.y
JetReconstruction.pz(recoparticle::ReconstructedParticle) = recoparticle.momentum.z
JetReconstruction.energy(recoparticle::ReconstructedParticle) = recoparticle.energy

function JetReconstruction.EEjet(recoparticle::ReconstructedParticle)
    EEjet(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
          JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle))
end

function JetReconstruction.PseudoJet(recoparticle::ReconstructedParticle)
    PseudoJet(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
              JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle))
end

end
