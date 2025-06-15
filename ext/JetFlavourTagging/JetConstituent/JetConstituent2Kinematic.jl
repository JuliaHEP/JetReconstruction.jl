using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### Basic Kinematic (11)

# get_pt - Transverse momentum
# get_p - Total momentum
# get_e - Energy
# get_mass - Mass
# get_type - Particle type
# get_charge - Electric charge
# get_theta - Polar angle
# get_phi - Azimuthal angle
# get_y - Rapidity
# get_eta - Pseudorapidity
# get_Bz - Magnetic field component

const JetConstituents = StructVector{ReconstructedParticle, <:Any}
const JetConstituentsData = Vector{Float32}

"""
    get_pt(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the transverse momentum of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)

# Returns
A vector of vectors of transverse momentum values (sqrt(px^2 + py^2))
"""
function get_pt(jcs::Vector{<:JetConstituents})
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        Float32[@inbounds sqrt(mom_x[i]^2 + mom_y[i]^2) 
                for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_p(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the momentum magnitude of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of momentum magnitudes (sqrt(px^2 + py^2 + pz^2))
"""
function get_p(jcs::Vector{<:JetConstituents})
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        mom_z = jet_constituents.momentum.z
        Float32[@inbounds sqrt(mom_x[i]^2 + mom_y[i]^2 + mom_z[i]^2) 
                for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_e(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the energy of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of energy values
"""
function get_e(jcs::Vector{JetConstituents})
    return [jet_constituents.energy for jet_constituents in jcs]
end

"""
    get_type(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the PDG type of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of particle types (PDG codes/Particle IDs)
"""
function get_type(jcs::Vector{<:JetConstituents})
    return [jet_constituents.type for jet_constituents in jcs]
end

"""
    get_mass(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the mass of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of mass values
"""
function get_mass(jcs::Vector{<:JetConstituents})
    return [jet_constituents.mass for jet_constituents in jcs]
end

"""
    get_charge(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the charge of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of charge values
"""
function get_charge(jcs::Vector{<:JetConstituents})
    return [jet_constituents.charge for jet_constituents in jcs]
end

"""
    get_theta(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the polar angle of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of polar angle values
"""
function get_theta(jcs::Vector{<:JetConstituents})
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        mom_z = jet_constituents.momentum.z
        Float32[@inbounds(
            let x = mom_x[i], y = mom_y[i], z = mom_z[i]
                (x == 0.0f0 && y == 0.0f0 && z == 0.0f0) ? 0.0f0 : atan(sqrt(x^2 + y^2), z)
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_phi(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the azimuthal angle of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of azimuthal angle values
"""
function get_phi(jcs::Vector{<:JetConstituents})
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        Float32[@inbounds(
            let x = mom_x[i], y = mom_y[i]
                (x == 0.0f0 && y == 0.0f0) ? 0.0f0 : atan(y, x)
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_y(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}
    
Get the rapidity of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of rapidity values
"""
function get_y(jcs::Vector{<:JetConstituents})
    return [begin
        energies = jet_constituents.energy
        mom_z = jet_constituents.momentum.z
        Float32[@inbounds(
            let e = energies[i], pz = mom_z[i]
                0.5f0 * log((e + pz) / (e - pz))
            end
        ) for i in eachindex(energies)]
    end for jet_constituents in jcs]
end

"""
    get_eta(jcs::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the pseudorapidity of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of pseudorapidity values (eta = -ln(tan(theta/2)))
"""
function get_eta(jcs::Vector{<:JetConstituents})
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        mom_z = jet_constituents.momentum.z
        Float32[@inbounds(
            let x = mom_x[i], y = mom_y[i], z = mom_z[i]
                p = sqrt(x^2 + y^2 + z^2)
                if p == 0.0f0
                    0.0f0
                elseif p == abs(z)  # particle along beam axis
                    sign(z) * Inf32
                else
                    0.5f0 * log((p + z) / (p - z))
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_Bz(jcs::Vector{<:JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Calculate the magnetic field Bz for each particle based on track curvature and momentum.

# Arguments
- jcs: Vector of jet constituents
- tracks: Vector of track states (used to get the omega value)

# Returns 
A vector of vectors of Bz values.
"""
function get_Bz(jcs::Vector{<:JetConstituents},
                tracks::StructVector{EDM4hep.TrackState})
    
    a = 2.99792458e8 * 1e3 * 1e-15
    n_tracks = length(tracks)
    
    # If tracks is a StructVector, we can access omega column directly
    omega_values = tracks.omega
    
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        charges = jet_constituents.charge
        track_indices = jet_constituents.tracks
        
        Float32[@inbounds(
            let track_idx = track_indices[i].first
                if track_idx < n_tracks
                    pt = sqrt(mom_x[i]^2 + mom_y[i]^2)
                    omega_values[track_idx + 1] / a * pt * copysign(1.0f0, charges[i])
                else
                    -9.0f0
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end