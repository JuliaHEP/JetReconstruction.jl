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
                            jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        jet_csts = Float32[]
        
        jet = jets[i]
        px, py, pz = jet.px, jet.py, jet.pz
        theta_jet = (px == 0.0f0 && py == 0.0f0 && pz == 0.0f0) ? 0.0f0 : atan(sqrt(px^2 + py^2), pz)
        phi_jet = (px == 0.0f0 && py == 0.0f0) ? 0.0f0 : atan(py, px)

        csts = jcs[i]
        for jc in csts
            p_const_x = jc.momentum.x
            p_const_y = jc.momentum.y
            p_const_z = jc.momentum.z

            # rotate z by -phi_jet, then rotate y by -theta_jet
            # First rotate around z-axis by -phi_jet
            p_rot_x = p_const_x * cos(-phi_jet) - p_const_y * sin(-phi_jet)
            p_rot_y = p_const_x * sin(-phi_jet) + p_const_y * cos(-phi_jet)
            p_rot_z = p_const_z
            # Then rotate around y-axis by -theta_jet
            p_rot2_x = p_rot_x * cos(-theta_jet) - p_rot_z * sin(-theta_jet)
            p_rot2_z = p_rot_x * sin(-theta_jet) + p_rot_z * cos(-theta_jet)
            p_rot2_y = p_rot_y
            # Calculate theta in rotated frame
            theta_rel = (p_rot2_x == 0.0f0 && p_rot2_y == 0.0f0 && p_rot2_z == 0.0f0) ? 0.0f0 : atan(sqrt(p_rot2_x^2 + p_rot2_y^2), p_rot2_z)
            push!(jet_csts, theta_rel)
        end
        push!(result, jet_csts)
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
                            jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        jet_csts = Float32[]

        jet = jets[i]
        px, py, pz= jet.px, jet.py, jet.pz
        theta_jet = (px == 0.0f0 && py == 0.0f0 && pz == 0.0f0) ? 0.0f0 : atan(sqrt(px^2 + py^2), pz)
        phi_jet = (px == 0.0f0 && py == 0.0f0) ? 0.0f0 : atan(py, px)
        
        csts = jcs[i]
        for jc in csts
            p_const_x = jc.momentum.x
            p_const_y = jc.momentum.y
            p_const_z = jc.momentum.z

            # rotate z by -phi_jet, then rotate y by -theta_jet
            # First rotate around z-axis by -phi_jet
            p_rot_x = p_const_x * cos(-phi_jet) - p_const_y * sin(-phi_jet)
            p_rot_y = p_const_x * sin(-phi_jet) + p_const_y * cos(-phi_jet)
            p_rot_z = p_const_z
            # Then rotate around y-axis by -theta_jet
            p_rot2_x = p_rot_x * cos(-theta_jet) - p_rot_z * sin(-theta_jet)
            p_rot2_z = p_rot_x * sin(-theta_jet) + p_rot_z * cos(-theta_jet) # Not needed.
            p_rot2_y = p_rot_y
            # Calculate phi in rotated frame
            phi_rel = (p_rot2_x == 0.0f0 && p_rot2_y == 0.0f0) ? 0.0f0 : atan(p_rot2_y, p_rot2_x)
            push!(jet_csts, phi_rel)
        end
        
        push!(result, jet_csts)
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
function get_thetaphirel_cluster(jets::Vector{EEJet}, 
                            jcs::Vector{JetConstituents})
    result = Tuple{Vector{JetConstituentsData}, Vector{JetConstituentsData}}(Vector{JetConstituentsData}(), Vector{JetConstituentsData}())
    
    for i in eachindex(jets)
        theta_csts = Float32[]
        phi_csts = Float32[]

        jet = jets[i]
        px, py, pz= jet.px, jet.py, jet.pz
        theta_jet = (px == 0.0f0 && py == 0.0f0 && pz == 0.0f0) ? 0.0f0 : atan(sqrt(px^2 + py^2), pz)
        phi_jet = (px == 0.0f0 && py == 0.0f0) ? 0.0f0 : atan(py, px)
        
        csts = jcs[i]
        for jc in csts
            p_const_x = jc.momentum.x
            p_const_y = jc.momentum.y
            p_const_z = jc.momentum.z

            # rotate z by -phi_jet, then rotate y by -theta_jet
            # First rotate around z-axis by -phi_jet
            p_rot_x = p_const_x * cos(-phi_jet) - p_const_y * sin(-phi_jet)
            p_rot_y = p_const_x * sin(-phi_jet) + p_const_y * cos(-phi_jet)
            p_rot_z = p_const_z
            # Then rotate around y-axis by -theta_jet
            p_rot2_x = p_rot_x * cos(-theta_jet) - p_rot_z * sin(-theta_jet)
            p_rot2_z = p_rot_x * sin(-theta_jet) + p_rot_z * cos(-theta_jet) # Not needed.
            p_rot2_y = p_rot_y
            # Calculate in rotated frame
            theta_rel = (p_rot2_x == 0.0f0 && p_rot2_y == 0.0f0 && p_rot2_z == 0.0f0) ? 0.0f0 : atan(sqrt(p_rot2_x^2 + p_rot2_y^2), p_rot2_z)
            phi_rel = (p_rot2_x == 0.0f0 && p_rot2_y == 0.0f0) ? 0.0f0 : atan(p_rot2_y, p_rot2_x)
            push!(theta_csts, theta_rel)
            push!(phi_csts, phi_rel)
        end

        push!(result[1], theta_csts)
        push!(result[2], phi_csts)
    end
    
    return result
end