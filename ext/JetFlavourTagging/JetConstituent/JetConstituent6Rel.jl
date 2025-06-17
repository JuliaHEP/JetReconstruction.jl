using EDM4hep
using JetReconstruction
using LoopVectorization
using StructArrays: StructVector

### Relative Kinematics (4)

# get_erel_cluster - Relative energy for clustered jets
# get_erel_log_cluster - Log of relative energy for clustered jets
# get_thetarel_cluster - Relative polar angle for clustered jets
# get_phirel_cluster - Relative azimuthal angle for clustered jets
# get_theta_phi_rel_cluster - Combined relative angles for clustered jets as they use the same logic


"""
    get_erel_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative energy (E_const/E_jet) for each constituent particle in clustered jets.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets.
- `jcs::Vector{<:JetConstituents}`: Vector of jet constituents corresponding to the jets.

# Returns
Vector containing relative energy values for each constituent in the jets.
"""
function get_erel_log_cluster(jets::Vector{EEJet}, 
                              jcs::Vector{<:JetConstituents})
    n_jets = length(jets)
    res = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        e_jet = jets[i].E
        csts = jcs[i]
        energies = csts.energy
        n_constituents = length(energies)
        jet_csts = Vector{Float32}(undef, n_constituents)
        
        if e_jet > 0.0f0
            inv_e_jet = 1.0f0 / e_jet
            @turbo for j in 1:n_constituents
                jet_csts[j] = energies[j] * inv_e_jet
            end
        else
            @turbo for j in 1:n_constituents
                jet_csts[j] = 0.0f0
            end
        end
        
        res[i] = jet_csts
    end
    
    return res
end

"""
    get_erel_log_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate log of relative energy (log(E_const/E_jet)) for each constituent particle in clustered jets.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jcs::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Vector containing log of relative energy values for each constituent in the jets
"""
function get_erel_log_cluster(jets::Vector{EEJet}, 
                              jcs::Vector{<:JetConstituents})
    n_jets = length(jets)
    res = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        e_jet = jets[i].E
        csts = jcs[i]
        energies = csts.energy
        n_constituents = length(energies)
        jet_csts = Vector{Float32}(undef, n_constituents)
        
        if e_jet > 0.0f0
            inv_e_jet = 1.0f0 / e_jet
            @turbo for j in 1:n_constituents
                jet_csts[j] = log10(energies[j] * inv_e_jet)
            end
        else
            @turbo for j in 1:n_constituents
                jet_csts[j] = 0.0f0
            end
        end
        
        res[i] = jet_csts
    end
    
    return res
end

"""
    get_thetarel_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative theta angle between constituent particle and clustered jet axis.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jcs::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Vector containing relative theta angle values for each constituent in the jets
"""
function get_thetarel_cluster(jets::Vector{EEJet}, 
                              jcs::Vector{<:JetConstituents})
    n_jets = length(jets)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        jet = jets[i]
        px, py, pz = jet.px, jet.py, jet.pz
        
        # Pre-compute jet angles
        pt_jet = sqrt(px^2 + py^2)
        theta_jet = atan(pt_jet, pz)
        phi_jet = atan(py, px)
        
        # Pre-compute trig values using sincos
        sin_phi, cos_phi = sincos(-phi_jet)
        sin_theta, cos_theta = sincos(-theta_jet)
        
        csts = jcs[i]
        mom_x = csts.momentum.x
        mom_y = csts.momentum.y
        mom_z = csts.momentum.z
        n_constituents = length(mom_x)
        jet_csts = Vector{Float32}(undef, n_constituents)
        
        @turbo for j in 1:n_constituents
            # First rotation
            p_rot_x = mom_x[j] * cos_phi - mom_y[j] * sin_phi
            p_rot_y = mom_x[j] * sin_phi + mom_y[j] * cos_phi
            
            # Second rotation
            p_rot2_x = p_rot_x * cos_theta - mom_z[j] * sin_theta
            p_rot2_z = p_rot_x * sin_theta + mom_z[j] * cos_theta

            pt_rot_sq = p_rot2_x^2 + p_rot_y^2
            pt_rot = sqrt(pt_rot_sq)
            jet_csts[j] = atan(pt_rot, p_rot2_z)
        end
        
        result[i] = jet_csts
    end
    
    return result
end

"""
    get_phirel_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative phi angle between constituent particle and clustered jet axis.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jcs::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Vector containing relative phi angle values for each constituent in the jets
"""
function get_phirel_cluster(jets::Vector{EEJet}, 
                            jcs::Vector{<:JetConstituents})
    n_jets = length(jets)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        jet = jets[i]
        px, py, pz = jet.px, jet.py, jet.pz
        
        # Pre-compute jet angles
        pt_jet = sqrt(px^2 + py^2)
        theta_jet = atan(pt_jet, pz)
        phi_jet = atan(py, px)
        
        # Pre-compute trig values using sincos
        sin_phi, cos_phi = sincos(-phi_jet)
        sin_theta, cos_theta = sincos(-theta_jet)
        
        csts = jcs[i]
        mom_x = csts.momentum.x
        mom_y = csts.momentum.y
        mom_z = csts.momentum.z
        n_constituents = length(mom_x)
        jet_csts = Vector{Float32}(undef, n_constituents)
        
        @turbo for j in 1:n_constituents
            # First rotation around z-axis by -phi_jet
            p_rot_x = mom_x[j] * cos_phi - mom_y[j] * sin_phi
            p_rot_y = mom_x[j] * sin_phi + mom_y[j] * cos_phi
            
            # Second rotation around y-axis by -theta_jet
            p_rot2_x = p_rot_x * cos_theta - mom_z[j] * sin_theta
            
            # Calculate phi in rotated frame
            jet_csts[j] = atan(p_rot_y, p_rot2_x)
        end
        
        result[i] = jet_csts
    end
    
    return result
end

"""
    get_thetaphirel_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative theta and phi angles between constituent particles and clustered jet axis.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jcs::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Tuple of Vectors containing relative theta and phi angle values for each constituent in the jets
"""
function get_thetarel_phirel_cluster(jets::Vector{EEJet}, 
                                     jcs::Vector{<:JetConstituents})
    n_jets = length(jets)
    theta_result = Vector{Vector{Float32}}(undef, n_jets)
    phi_result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        jet = jets[i]
        px, py, pz = jet.px, jet.py, jet.pz
        
        # Pre-compute jet angles
        pt_jet = sqrt(px^2 + py^2)
        theta_jet = atan(pt_jet, pz)
        phi_jet = atan(py, px)
        
        # Pre-compute trig values using sincos
        sin_phi, cos_phi = sincos(-phi_jet)
        sin_theta, cos_theta = sincos(-theta_jet)
        
        csts = jcs[i]
        mom_x = csts.momentum.x
        mom_y = csts.momentum.y
        mom_z = csts.momentum.z
        n_constituents = length(mom_x)
        
        jet_theta = Vector{Float32}(undef, n_constituents)
        jet_phi = Vector{Float32}(undef, n_constituents)
        
        @turbo for j in 1:n_constituents
            # First rotation around z-axis by -phi_jet
            p_rot_x = mom_x[j] * cos_phi - mom_y[j] * sin_phi
            p_rot_y = mom_x[j] * sin_phi + mom_y[j] * cos_phi
            
            # Second rotation around y-axis by -theta_jet
            p_rot2_x = p_rot_x * cos_theta - mom_z[j] * sin_theta
            p_rot2_z = p_rot_x * sin_theta + mom_z[j] * cos_theta
            p_rot2_y = p_rot_y
            
            # Calculate both theta and phi in rotated frame
            pt_rot = sqrt(p_rot2_x^2 + p_rot2_y^2)
            jet_theta[j] = atan(pt_rot, p_rot2_z)
            jet_phi[j] = atan(p_rot2_y, p_rot2_x)
        end
        
        theta_result[i] = jet_theta
        phi_result[i] = jet_phi
    end
    
    return (theta_result, phi_result)
end