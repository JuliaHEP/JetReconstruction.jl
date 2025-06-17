using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### Particle Identification (5)

# get_isMu - Check if constituent is muon
# get_isEl - Check if constituent is electron
# get_isChargedHad - Check if constituent is charged hadron
# get_isGamma - Check if constituent is photon
# get_isNeutralHad - Check if constituent is neutral hadron

"""
    get_isMu(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a muon.

# Arguments
- jcs: Vector of jet constituents 

# Returns
Vector of vectors of is muon boolean values as Float32.
"""
function get_isMu(jcs::Vector{<:JetConstituents})
    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        charges = jet_constituents.charge
        masses = jet_constituents.mass
        n_particles = length(charges)
        
        is_mu = Vector{Float32}(undef, n_particles)
        
        @simd for i in 1:n_particles
            charge_check = abs(charges[i]) > 0
            mass_check = abs(masses[i] - 0.105658f0) < 1.0f-3
            is_mu[i] = (charge_check & mass_check) ? 1.0f0 : 0.0f0
        end
        
        result[j] = is_mu
    end
    
    return result
end

"""
    get_isEl(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is an electron.

# Arguments
- jcs: Vector of jet constituents 

# Returns
Vector of vectors of is electron boolean values as Float32.
"""
function get_isEl(jcs::Vector{<:JetConstituents})

    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        charges = jet_constituents.charge
        masses = jet_constituents.mass
        n_particles = length(charges)
        
        is_mu = Vector{Float32}(undef, n_particles)
        
        @simd for i in 1:n_particles
            charge_check = abs(charges[i]) > 0
            mass_check = abs(masses[i] - 0.000510999f0) < 1.0f-5
            is_mu[i] = (charge_check & mass_check) ? 1.0f0 : 0.0f0
        end
        
        result[j] = is_mu
    end
    
    return result
end

"""
    get_isChargedHad(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a charged hadron.
    
# Arguments
- jcs: Vector of jet constituents 

# Returns
Vector of vectors of is charged hadron boolean values as Float32.
"""
function get_isChargedHad(jcs::Vector{JetConstituents})
    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        charges = jet_constituents.charge
        masses = jet_constituents.mass
        n_particles = length(charges)
        
        is_mu = Vector{Float32}(undef, n_particles)
        
        @simd for i in 1:n_particles
            charge_check = abs(charges[i]) > 0
            mass_check = abs(masses[i] - 0.13957f0) < 1.0f-3
            is_mu[i] = (charge_check & mass_check) ? 1.0f0 : 0.0f0
        end
        
        result[j] = is_mu
    end
    
    return result
end

"""
    get_isGamma(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a photon (gamma) (PDG 22).

# Arguments
- jcs: Vector of jet constituents 

# Returns
Vector of vectors of is photon boolean values as Float32.
"""
function get_isGamma(jcs::Vector{<:JetConstituents})
    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        types = jet_constituents.type
        n_particles = length(types)
        
        is_gamma = Vector{Float32}(undef, n_particles)
        
        @simd for i in 1:n_particles
            is_gamma[i] = types[i] == 22 ? 1.0f0 : 0.0f0
        end
        
        result[j] = is_gamma
    end
    
    return result
end

"""
    get_isNeutralHad(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a neutral hadron (PDG 130).

# Arguments
- jcs: Vector of jet constituents 

# Returns
Vector of vectors of is neutral hadron boolean values as Float32.
"""
function get_isNeutralHad(jcs::Vector{<:JetConstituents})
    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        types = jet_constituents.type
        n_particles = length(types)
        
        is_gamma = Vector{Float32}(undef, n_particles)
        
        @simd for i in 1:n_particles
            is_gamma[i] = types[i] == 130 ? 1.0f0 : 0.0f0
        end
        
        result[j] = is_gamma
    end
    
    return result
end