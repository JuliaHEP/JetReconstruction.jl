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

import Base.eltype
"""
    eltype(recoparticle::ReconstructedParticle)

Return the underlying numerical type of the reconstructed particle
"""
eltype(recoparticle::ReconstructedParticle) = typeof(recoparticle.energy)

"""
    JetReconstruction.EEJet{T}(recoparticle::ReconstructedParticle) where {T <: Real}

Construct an EEJet{T} from a ReconstructedParticle.
"""
function JetReconstruction.EEJet{T}(recoparticle::ReconstructedParticle;
                                    cluster_hist_index::Int = 0) where {T <: Real}
    EEJet{T}(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
             JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle);
             cluster_hist_index)
end

"""
    JetReconstruction.EEJet(recoparticle::ReconstructedParticle)

Construct an EEJet from a ReconstructedParticle, taking the numerical type from
the incoming recoparticle.
"""
function JetReconstruction.EEJet(recoparticle::ReconstructedParticle;
                                 cluster_hist_index::Int = 0)
    T = eltype(recoparticle)
    EEJet{T}(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
             JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle);
             cluster_hist_index)
end

"""
    JetReconstruction.PseudoJet{T}(recoparticle::ReconstructedParticle) where {T <: Real}

Construct an PseudoJet{T} from a ReconstructedParticle.
"""
function JetReconstruction.PseudoJet{T}(recoparticle::ReconstructedParticle;
                                        cluster_hist_index::Int = 0) where {T <: Real}
    PseudoJet{T}(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
                 JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle);
                 cluster_hist_index)
end

"""
    JetReconstruction.PseudoJet(recoparticle::ReconstructedParticle)

Construct an PseudoJet from a ReconstructedParticle, taking the numerical type from
the incoming recoparticle.
"""
function JetReconstruction.PseudoJet(recoparticle::ReconstructedParticle;
                                     cluster_hist_index::Int = 0)
    T = eltype(recoparticle)
    PseudoJet{T}(JetReconstruction.px(recoparticle), JetReconstruction.py(recoparticle),
                 JetReconstruction.pz(recoparticle), JetReconstruction.energy(recoparticle);
                 cluster_hist_index)
end

end
