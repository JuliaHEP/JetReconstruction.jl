## Jet class constructors from EMD4hep objects

module EDM4hepJets

using JetReconstruction
using EDM4hep

# Make a union so that we can consider MC and Reco particles to
# be the same types for the following methods
const EDM4hepParticle = Union{MCParticle, ReconstructedParticle}

"""
    JetReconstruction.px(particle::EDM4hepParticle)

Return the x component of the momentum of a EDM4hepParticle.
"""
JetReconstruction.px(particle::EDM4hepParticle) = particle.momentum.x

"""
    JetReconstruction.py(particle::EDM4hepParticle)

Return the y component of the momentum of a EDM4hepParticle.
"""
JetReconstruction.py(particle::EDM4hepParticle) = particle.momentum.y

"""
    JetReconstruction.pz(particle::EDM4hepParticle)

Return the z component of the momentum of a EDM4hepParticle.
"""
JetReconstruction.pz(particle::EDM4hepParticle) = particle.momentum.z

"""
    JetReconstruction.energy(particle::EDM4hepParticle)

Return the energy component of a EDM4hepParticle's four vector.
"""
JetReconstruction.energy(particle::EDM4hepParticle) = particle.energy

"""
    JetReconstruction.EEJet(particle::EDM4hepParticle)

Construct an EEJet from a EDM4hepParticle.
"""
function JetReconstruction.EEJet(particle::EDM4hepParticle;
                                 cluster_hist_index::Int = 0)
    EEJet(JetReconstruction.px(particle), JetReconstruction.py(particle),
          JetReconstruction.pz(particle), JetReconstruction.energy(particle);
          cluster_hist_index)
end

"""
    JetReconstruction.PseudoJet(particle::EDM4hepParticle)

Construct an PseudoJet from a EDM4hepParticle.
"""
function JetReconstruction.PseudoJet(particle::EDM4hepParticle;
                                     cluster_hist_index::Int = 0)
    PseudoJet(JetReconstruction.px(particle), JetReconstruction.py(particle),
              JetReconstruction.pz(particle), JetReconstruction.energy(particle);
              cluster_hist_index)
end

end
