using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### Impact Parameters and Jet Distance (12)

## 2D Impact Parameter
# get_Sip2dVal_clusterV - Vectorized version 2D signed impact parameter value for clustered jets
# get_Sip2dSig - 2D impact parameter significance

# 3D Impact Parameter
# get_Sip3dVal_clusterV - Vectorized version 3D signed impact parameter value for clustered jets
# get_Sip3dSig - 3D impact parameter significance

# Jet Distance
# get_JetDistVal_clusterV - Vectorized version jet distance value for clustered jets
# get_JetDistSig - Jet distance significance

"""
    get_Sip2dVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                        D0::Vector{JetConstituentsData},
                        phi0::Vector{JetConstituentsData},
                        Bz::Float32) -> Vector{JetConstituentsData}

Calculate the 2D signed impact parameter value for each particle relative to the jet axis.
This is a Julia implementation of the C++ function get_Sip2dVal_clusterV.

# Arguments
- jets: Vector of EEJet objects representing jets
- D0: Vector of vectors containing D0 values (transverse impact parameters)
- phi0: Vector of vectors containing phi0 values (azimuthal angles at impact point)
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of 2D signed impact parameter values (one vector per jet)
"""
function get_Sip2dVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                            D0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32)
    n_jets = length(jets)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        px = Float32(jets[i].px)
        py = Float32(jets[i].py)
        d0_vals = D0[i]
        phi_vals = phi0[i]
        n_constituents = length(d0_vals)
        
        sip2d_values = Vector{Float32}(undef, n_constituents)
        
        @turbo for j in 1:n_constituents
            d0_val = d0_vals[j]
            phi_val = phi_vals[j]

            sin_phi, cos_phi = sincos(phi_val)

            d0x = -d0_val * sin_phi
            d0y = d0_val * cos_phi
            
            dot_product = d0x * px + d0y * py
            
            abs_d0 = abs(d0_val)
            sign_dot = sign(dot_product)
            signed_ip = sign_dot * abs_d0
            
            is_valid = Float32(d0_val != -9.0f0)
            sip2d_values[j] = is_valid * signed_ip + (1.0f0 - is_valid) * (-9.0f0)
        end
        
        result[i] = sip2d_values
    end
    
    return result
end


"""
    get_btagSip2dVal(jets::Vector{JetReconstruction.EEJet},
                    pfcand_dxy::Vector{JetConstituentsData},
                    pfcand_phi0::Vector{JetConstituentsData},
                    Bz::Float32) -> Vector{JetConstituentsData}

Call the implementation function get_Sip2dVal_clusterV
"""
function get_btagSip2dVal(jets::Vector{JetReconstruction.EEJet},
    pfcand_dxy::Vector{JetConstituentsData},
    pfcand_phi0::Vector{JetConstituentsData},
    Bz::Float32)
    # Simply call the implementation function
    return get_Sip2dVal_clusterV(jets, pfcand_dxy, pfcand_phi0, Bz)
end

"""
    get_Sip2dSig(Sip2dVals::Vector{JetConstituentsData},
                err2_D0::Vector{JetConstituentsData}) -> Vector{JetConstituentsData}

Calculate the 2D signed impact parameter significance for each particle.
This is a Julia implementation of the C++ function get_Sip2dSig.

# Arguments
- Sip2dVals: Vector of vectors containing 2D signed impact parameter values
- err2_D0: Vector of vectors containing squared errors of the D0 values

# Returns
Vector of vectors of 2D signed impact parameter significances (one vector per jet)
"""
function get_Sip2dSig(Sip2dVals::Vector{JetConstituentsData},
                    err2_D0::Vector{JetConstituentsData})
    n_jets = length(Sip2dVals)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        n_constituents = length(Sip2dVals[i])
        sig_values = Vector{Float32}(undef, n_constituents)
        sip_vals = Sip2dVals[i]
        err_vals = err2_D0[i]
        
        @simd for j in 1:n_constituents
            err_val = err_vals[j]
            sip_val = sip_vals[j]
            
            valid = err_val > 0.0f0
            sqrt_err = sqrt(max(err_val, eps(Float32)))  # Avoid sqrt of negative
            sig = sip_val / sqrt_err
            
            sig_values[j] = valid ? sig : -9.0f0
        end
        
        result[i] = sig_values
    end
    
    return result
end

"""
    get_btagSip2dSig(pfcand_btagSip2dVal::Vector{JetConstituentsData},
                    pfcand_dxydxy::Vector{JetConstituentsData}) -> Vector{JetConstituentsData}

Call the implementation function get_Sip2dSig
"""
function get_btagSip2dSig(pfcand_btagSip2dVal::Vector{JetConstituentsData},
                        pfcand_dxydxy::Vector{JetConstituentsData})
    # Simply call the implementation function
    return get_Sip2dSig(pfcand_btagSip2dVal, pfcand_dxydxy)
end

"""
    get_Sip3dVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                        D0::Vector{JetConstituentsData},
                        Z0::Vector{JetConstituentsData},
                        phi0::Vector{JetConstituentsData},
                        Bz::Float32) -> Vector{JetConstituentsData}

Calculate the 3D signed impact parameter value for each particle relative to the jet axis.
"""
function get_Sip3dVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32)
    n_jets = length(jets)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        px = Float32(jets[i].px)
        py = Float32(jets[i].py)
        pz = Float32(jets[i].pz)
        d0_vals = D0[i]
        z0_vals = Z0[i]
        phi_vals = phi0[i]
        n_constituents = length(d0_vals)

        cprojs = Vector{Float32}(undef, n_constituents)
        
        @turbo for j in 1:n_constituents
            d0_val = d0_vals[j]
            z0_val = z0_vals[j]
            phi_val = phi_vals[j]

            sin_phi = sin(phi_val)
            cos_phi = cos(phi_val)

            dx = -d0_val * sin_phi
            dy = d0_val * cos_phi
            dz = z0_val
            
            dot_prod = dx * px + dy * py + dz * pz
            
            magnitude = sqrt(d0_val * d0_val + z0_val * z0_val)
            sign_dot = sign(dot_prod)
            signed_ip = sign_dot * magnitude
            
            is_valid = Float32(d0_val != -9.0f0)
            cprojs[j] = is_valid * signed_ip + (1.0f0 - is_valid) * (-9.0f0)
        end
        
        result[i] = cprojs
    end
    
    return result
end

function get_btagSip3dVal(jets::Vector{JetReconstruction.EEJet},
    pfcand_dxy::Vector{JetConstituentsData},
    pfcand_dz::Vector{JetConstituentsData},
    pfcand_phi0::Vector{JetConstituentsData},
    Bz::Float32)
# Simply call the implementation function
    return get_Sip3dVal_clusterV(jets, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)
end

"""
    get_Sip3dSig(Sip3dVals::Vector{JetConstituentsData},
                err2_D0::Vector{JetConstituentsData},
                err2_Z0::Vector{JetConstituentsData}) -> Vector{JetConstituentsData}

Calculate the 3D signed impact parameter significance (value/error) for each particle.
"""
function get_Sip3dSig(Sip3dVals::Vector{JetConstituentsData},
                     err2_D0::Vector{JetConstituentsData},
                     err2_Z0::Vector{JetConstituentsData})
    n_jets = length(Sip3dVals)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        n_constituents = length(Sip3dVals[i])
        sigs = Vector{Float32}(undef, n_constituents)
        sip_vals = Sip3dVals[i]
        err_d0 = err2_D0[i]
        err_z0 = err2_Z0[i]
        
        @simd for j in 1:n_constituents
            err_d0_val = err_d0[j]
            err_z0_val = err_z0[j]
            sip_val = sip_vals[j]
            
            # Branchless computation
            valid = err_d0_val > 0.0f0
            err_sum = err_d0_val + err_z0_val
            sqrt_err = sqrt(max(err_sum, eps(Float32)))
            sig = sip_val / sqrt_err
            
            sigs[j] = valid ? sig : -9.0f0
        end
        
        result[i] = sigs
    end
    
    return result
end

function get_btagSip3dSig(Sip3dVals::Vector{JetConstituentsData},
                        err2_D0::Vector{JetConstituentsData},
                        err2_Z0::Vector{JetConstituentsData})
    # Simply call the implementation function
    return get_Sip3dSig(Sip3dVals, err2_D0, err2_Z0)
end

"""
    get_JetDistVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                            jcs::Vector{<:JetConstituents},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the jet distance value for each particle, measuring the distance between
the point of closest approach and the jet axis.
"""
function get_JetDistVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                                jcs::Vector{<:JetConstituents},
                                D0::Vector{JetConstituentsData},
                                Z0::Vector{JetConstituentsData},
                                phi0::Vector{JetConstituentsData},
                                Bz::Float32)
    n_jets = length(jets)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    for i in 1:n_jets
        px_jet, py_jet, pz_jet = jets[i].px, jets[i].py, jets[i].pz
        ct = jcs[i]
        n_constituents = length(ct)
        tmp = Vector{Float32}(undef, n_constituents)
        
        for j in 1:n_constituents
            if D0[i][j] != -9.0f0
                d0_val = D0[i][j]
                z0_val = Z0[i][j]
                
                # Use sincos for efficiency
                sin_phi, cos_phi = sincos(phi0[i][j])
                
                # Impact parameter vector
                dx = -d0_val * sin_phi
                dy = d0_val * cos_phi
                dz = z0_val
                
                # Constituent momentum
                px_ct = ct[j].momentum.x
                py_ct = ct[j].momentum.y
                pz_ct = ct[j].momentum.z
                
                # Cross product: n = p_ct × p_jet
                nx = py_ct * pz_jet - pz_ct * py_jet
                ny = pz_ct * px_jet - px_ct * pz_jet
                nz = px_ct * py_jet - py_ct * px_jet
                
                # Normalize
                n_mag = sqrt(nx^2 + ny^2 + nz^2)
                inv_n_mag = 1.0f0 / max(n_mag, eps(Float32))
                nx *= inv_n_mag
                ny *= inv_n_mag
                nz *= inv_n_mag
                
                # Distance (r_jet = [0,0,0], so we just need n·d)
                tmp[j] = nx * dx + ny * dy + nz * dz
            else
                tmp[j] = -9.0f0
            end
        end
        
        result[i] = tmp
    end
    
    return result
end

function get_btagJetDistVal(jets::Vector{JetReconstruction.EEJet},
                            jcs::Vector{<:JetConstituents},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32)
    # Simply call the implementation function
    return get_JetDistVal_clusterV(jets, jcs, D0, Z0, phi0, Bz)
end


"""
    get_JetDistSig(JetDistVal::Vector{JetConstituentsData},
                    err2_D0::Vector{JetConstituentsData},
                    err2_Z0::Vector{JetConstituentsData}) -> Vector{JetConstituentsData}

Calculate the jet distance significance (value/error) for each particle.
"""
function get_JetDistSig(JetDistVal::Vector{JetConstituentsData},
                        err2_D0::Vector{JetConstituentsData},
                        err2_Z0::Vector{JetConstituentsData})
    n_jets = length(JetDistVal)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for i in 1:n_jets
        n_constituents = length(JetDistVal[i])
        tmp = Vector{Float32}(undef, n_constituents)
        jet_vals = JetDistVal[i]
        err_d0 = err2_D0[i]
        err_z0 = err2_Z0[i]
        
        @simd for j in 1:n_constituents
            err_d0_val = err_d0[j]
            err_z0_val = err_z0[j]
            jet_val = jet_vals[j]
            
            # Branchless computation
            valid = err_d0_val > 0.0f0
            err_sum = err_d0_val + err_z0_val
            sqrt_err = sqrt(max(err_sum, eps(Float32)))
            sig = jet_val / sqrt_err
            
            tmp[j] = valid ? sig : -9.0f0
        end
        
        result[i] = tmp
    end
    
    return result
end

function get_btagJetDistSig(JetDistVal::Vector{JetConstituentsData},
                            err2_D0::Vector{JetConstituentsData},
                            err2_Z0::Vector{JetConstituentsData})
    # Simply call the implementation function
    return get_JetDistSig(JetDistVal, err2_D0, err2_Z0)
end
