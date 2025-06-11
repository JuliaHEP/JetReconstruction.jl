using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### TODO's:
# 1) Fix the covariance matrix extraction to avoid out-of-bounds errors. (DONE, Jun 9, 2025)
# 2) Move class type alias in line 11-12 to JetFlavourTagging file. But requires more checks.
# 3) For Vector structs and operations, maybe those file can go to another vector helper module, or I can use julia buildin Vector type.
# 4) Fix the kinematic variable sections (Currently working on, expect Jun 11, 2025)
# 5) Fix the vertex extraction from MCParticles 
#   a) prepare the vertex in the main file and pass it to here (currently working on, expect Jun 13, 2025)
#   b) MTOF, dNdx (done)
#   c) dNdx based on EDM4hep version 
#   d) kinematic variables (currently working on, expect Jun 16, 2025)
# 6) Optimize the b-tagging section. (expected Jun 17, 2025)
# 7) Unit tests! (Extract one event!)

# Define type aliases for clarity 
const JetConstituents = StructVector{ReconstructedParticle}
const JetConstituentsData = Vector{Float32}

#############################################################################
# Vector structs and operations
#############################################################################

"""
    Vector2(x::Float64, y::Float64) -> Vector2

Create a 2D vector with x and y components.

Fields: 
- `X`: The x component of the vector.
- `Y`: The y component of the vector.
"""
struct Vector2
    X::Float64
    Y::Float64
end

"""
    Vector3(x::Float64, y::Float64, z::Float64) -> Vector3

Create a 3D vector with x, y, and z components.

Fields:
- `X`: The x component of the vector.
- `Y`: The y component of the vector.
- `Z`: The z component of the vector.
"""
struct Vector3
    X::Float64
    Y::Float64
    Z::Float64
end

function v3_cross(v::Vector3, w::Vector3)
    return Vector3(
        v.Y*w.Z - v.Z*w.Y,
        v.Z*w.X - v.X*w.Z,
        v.X*w.Y - v.Y*w.X,
    )
end

"""
    v3_dot(v::Vector3, w::Vector3) -> Float64

Calculate the dot product of two 3D vectors.

# Returns
The dot product of the two vectors.
"""
function v3_dot(v::Vector3, w::Vector3)
    return v.X*w.X + v.Y*w.Y + v.Z*w.Z
end

"""
    v3_norm(v::Vector3) -> Float64

Calculate the norm (magnitude) of a 3D vector.

# Returns
The magnitude of the vector.
"""
function v3_norm(v::Vector3)
    return sqrt(v.X^2 + v.Y^2 + v.Z^2)
end

"""
    v3_unit(v::Vector3) -> Vector3

Calculate the unit vector of a 3D vector.

# Returns
The unit vector in the same direction as the input vector.
"""
function v3_unit(v::Vector3)
    n = v3_norm(v)
    return Vector3(v.X/n, v.Y/n, v.Z/n)
end

#################################################################################
# Kinematic variables
#################################################################################
"""
    get_Bz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Calculate the magnetic field Bz for each particle based on track curvature and momentum.
# Returns 
A vector of vectors of Bz values (one vector per jet, one value per constituent).
"""
function get_Bz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState})
    # Constants
    c_light = 2.99792458e8  # speed of light in m/s
    a = c_light * 1e3 * 1e-15  # conversion factor for omega [1/mm]
    
    result = Vector{JetConstituentsData}()
    
    for constituents in jcs
        bz_values = JetConstituentsData()
        
        for p in constituents
            # Check if particle has associated tracks through the relation
            if isdefined(p, :tracks) && !isnothing(p.tracks) && !isempty(p.tracks)
                # Get the first track (most relevant for Bz calculation)
                track_idx = p.tracks[1]
                
                if track_idx <= length(tracks)
                    track = tracks[track_idx]
                    pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
                    
                    # Calculate Bz based on track curvature (omega) and momentum
                    charge_sign = p.charge > 0.0 ? 1.0 : -1.0
                    bz_val = track.omega / a * pt * charge_sign
                    push!(bz_values, bz_val)
                else
                    push!(bz_values, -9.0f0)
                end
            else
                push!(bz_values, -9.0f0)
            end
        end
        
        push!(result, bz_values)
    end
    
    return result
end

"""
    get_mass(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the mass of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)

# Returns
A vector of vectors of mass values
"""
function get_mass(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        mass_values = JetConstituentsData()
        for p in jet_constituents
            mass = p.mass
            push!(mass_values, mass)
        end
        push!(result, mass_values)
    end
    return result
end

"""
    get_eta(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the pseudorapidity of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)

# Returns
A vector of vectors of pseudorapidity values (eta = -ln(tan(theta/2)))
"""
function get_eta(jcs::Vector{JetConstituents})
    theta_result = get_theta(jcs)
    result = -ln.(tan.(theta_result ./ 2))
    return result
end

"""
    get_pt(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the transverse momentum of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)

# Returns
A vector of vectors of transverse momentum values (sqrt(px^2 + py^2))
"""
function get_pt(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        pt_values = JetConstituentsData()
        for p in jet_constituents
            pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
            push!(pt_values, pt)
        end
        push!(result, pt_values)
    end
    return result
end

"""
    get_p(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the momentum magnitude of each particle in each jet.

# Arguments
- jcs: Vector of jet constituents

# Returns
A vector of vectors of momentum magnitudes (sqrt(px^2 + py^2 + pz^2))
"""
function get_p(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        p_values = JetConstituentsData()
        for p in jet_constituents
            momentum = sqrt(p.momentum.x^2 + p.momentum.y^2 + p.momentum.z^2)
            push!(p_values, momentum)
        end
        push!(result, p_values)
    end
    return result
end

"""
    get_e(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the energy of each particle in each jet.
"""
function get_e(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        e_values = JetConstituentsData()
        for p in jet_constituents
            push!(e_values, p.energy)
        end
        push!(result, e_values)
    end
    return result
end

"""
    get_theta(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the polar angle of each particle in each jet.
"""
function get_theta(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        theta_values = JetConstituentsData()
        for p in jet_constituents
            theta = (p.momentum.x == 0.0 && p.momentum.y == 0.0 && p.momentum.z == 0.0) ? 0.0 : atan(sqrt(p.momentum.x^2 + p.momentum.y^2), p.momentum.z)
            push!(theta_values, theta)
        end
        push!(result, theta_values)
    end
    return result
end

"""
    get_phi(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the azimuthal angle of each particle in each jet.
"""
function get_phi(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        phi_values = JetConstituentsData()
        for p in jet_constituents
            phi = (p.momentum.x == 0.0 && p.momentum.y == 0.0) ? 0.0 : atan(p.momentum.y, p.momentum.x)
            push!(phi_values, phi)
        end
        push!(result, phi_values)
    end
    return result
end

"""
    get_y(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}
    
Get the rapidity of each particle in each jet.
"""
function get_y(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        y_values = JetConstituentsData()
        for p in jet_constituents
            y = 0.5 * log((p.energy + p.momentum.z) / (p.energy - p.momentum.z))
            push!(y_values, y)
        end
        push!(result, y_values)
    end
    return result
end

"""
    get_charge(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the charge of each particle in each jet.
"""
function get_charge(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        charge_values = JetConstituentsData()
        for p in jet_constituents
            push!(charge_values, p.charge)
        end
        push!(result, charge_values)
    end
    return result
end

"""
    get_type(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the PDG type of each particle in each jet.

# Returns
A vector of vectors of particle types (PDG codes/Particle IDs)
"""
function get_type(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    for jet_constituents in jcs
        type_values = JetConstituentsData()
        for p in jet_constituents
            push!(type_values, Float32(p.type))
        end
        push!(result, type_values)
    end
    return result
end

"""
    get_phi0(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the phi angle at the point of closest approach for each particle relative to vertex V.
This is a Julia implementation of the C++ function XPtoPar_phi.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of phi values (one vector per jet)
"""
function get_phi0(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                V::LorentzVector, Bz::Float32)

    cSpeed = 2.99792458e8 * 1.0e-9 
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)
    
    for jc in jcs
        phi_values = JetConstituentsData()
        for p in jc
            if p.tracks.first < tracks_len
                track = tracks[p.tracks.first + 1]
                
                D0_wrt0 = track.D0
                Z0_wrt0 = track.Z0
                phi0_wrt0 = track.phi

                X = [-D0_wrt0 * sin(phi0_wrt0), 
                    D0_wrt0 * cos(phi0_wrt0), 
                    Z0_wrt0]
                
                x = X .- [V.x, V.y, V.z]
                p3 = [p.momentum.x, p.momentum.y, p.momentum.z]
                a = -p.charge * Bz * cSpeed
                pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
                r2 = x[1]^2 + x[2]^2
                cross = x[1] * p3[2] - x[2] * p3[1]
                T = sqrt(pt^2 - 2 * a * cross + a^2 * r2)
                phi0 = atan((p3[2] - a * x[1]) / T, (p3[1] + a * x[2]) / T)

                push!(phi_values, Float32(phi0))
            else
                push!(phi_values, -9.0f0)
            end
        end
        push!(result, phi_values)
    end

    return result
end

"""
    get_dxy(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the transverse impact parameter dxy for each particle in each jet relative to vertex V.
Reference: FCCAnalyses c++ function XPtoPar_dxy, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of dxy values (one vector per jet)
"""
function get_dxy(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                V::LorentzVector, Bz::Float32)

    cSpeed = 2.99792458e8 * 1.0e-9 
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)
    
    for jc in jcs
        dxy_values = JetConstituentsData()
        for p in jc
            if p.tracks.first < tracks_len
                track = tracks[p.tracks.first + 1]
                
                D0_wrt0 = track.D0
                Z0_wrt0 = track.Z0
                phi0_wrt0 = track.phi

                X = [-D0_wrt0 * sin(phi0_wrt0), 
                    D0_wrt0 * cos(phi0_wrt0), 
                    Z0_wrt0]
                
                x = X .- [V.x, V.y, V.z]
                p3 = [p.momentum.x, p.momentum.y, p.momentum.z]
                a = -p.charge * Bz * cSpeed
                pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
                r2 = x[1]^2 + x[2]^2
                cross = x[1] * p3[2] - x[2] * p3[1]
                D = -9.0f0

                if (pt^2 - 2 * a * cross + a^2 * r2) > 0
                    T = sqrt(pt^2 - 2 * a * cross + a^2 * r2)
                    if pt < 10.0
                        D = (T - pt) / a
                    else
                        D = (-2 * cross + a * r2) / (T + pt)
                    end
                end
                push!(dxy_values, D)
            else
                push!(dxy_values, -9.0f0)
            end
        end
        push!(result, dxy_values)
    end
    return result
end

"""
    get_dz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the longitudinal impact parameter dz for each particle in each jet relative to vertex V.
Reference: FCCAnalyses c++ function XPtoPar_dz, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of dz values (one vector per jet)
"""
function get_dz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                V::LorentzVector, Bz::Float32)

    cSpeed = 2.99792458e8 * 1.0e-9 
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)
    
    for jc in jcs
        dz_values = JetConstituentsData()
        for p in jc
            if p.tracks.first < tracks_len
                track = tracks[p.tracks.first + 1]
                
                D0_wrt0 = track.D0
                Z0_wrt0 = track.Z0
                phi0_wrt0 = track.phi

                X = [-D0_wrt0 * sin(phi0_wrt0), 
                    D0_wrt0 * cos(phi0_wrt0), 
                    Z0_wrt0]
                
                x = X .- [V.x, V.y, V.z]
                p3 = [p.momentum.x, p.momentum.y, p.momentum.z]
                a = -p.charge * Bz * cSpeed
                pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
                C = a / (2 * pt)
                r2 = x[1]^2 + x[2]^2
                cross = x[1] * p3[2] - x[2] * p3[1]
                T = sqrt(pt^2 - 2 * a * cross + a^2 * r2)

                D = if pt < 10.0
                    (T - pt) / a
                else
                    (-2 * cross + a * r2) / (T + pt)
                end

                B = C * sqrt(max(r2 - D^2, 0.0) / (1 + 2 * C * D))
                if abs(B) > 1.0
                    B = sign(B)
                end

                st = asin(B) / C
                ct = p3[3] / pt

                dot = x[1] * p3[1] + x[2] * p3[2]
                z0 = if dot > 0.0
                    x[3] - ct * st
                else
                    x[3] + ct * st
                end

                push!(dz_values, z0)
            else
                push!(dz_values, -9.0f0)
            end
        end
        push!(result, dz_values)
    end
    return result
end

"""
    get_c(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the C for each particle in each jet.
Reference: FCCAnalyses c++ function XPtoPar_C, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of C values (one vector per jet)
"""
function get_c(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                Bz::Float32)

    cSpeed = 2.99792458e8 * 1.0e3 * 1.0e-15
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)
    
    for jc in jcs
        c_values = JetConstituentsData()
        for p in jc
            if p.tracks.first < tracks_len
                a = copysign(1.0, p.charge) * Bz * cSpeed
                pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
                C = a/(2 * pt)
                push!(c_values, C)
            else
                push!(c_values, -9.0f0)
            end
        end
        push!(result, c_values)
    end
    return result
end

"""
    get_ct(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the ct for each particle in each jet.
Reference: FCCAnalyses c++ function XPtoPar_C, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of ct values (one vector per jet)
"""
function get_ct(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                Bz::Float32)

    cSpeed = 2.99792458e8 * 1.0e-9
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)
    
    for jc in jcs
        ct_values = JetConstituentsData()
        for p in jc
            if p.tracks.first < tracks_len
                pt = sqrt(p.momentum.x^2 + p.momentum.y^2)
                ct = p.momentum.z / pt
                push!(ct_values, ct)
            else
                push!(ct_values, -9.0f0)
            end
        end
        push!(result, ct_values)
    end
    return result
end

#################################################################################
# Covariance matrix elements
#################################################################################4

"""
    get_dxydxy(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0 covariance (dxy/dxy) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxydxy values (one vector per jet)
"""
function get_dxydxy(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})

    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[1])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dphidxy(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-d0 covariance (dphi/dxy) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dphidxy values (one vector per jet)
"""
function get_dphidxy(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})

    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[2])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dphidphi(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi covariance (dphi/dphi) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_phi0_cov, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dphidphi values (one vector per jet)
"""
function get_dphidphi(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})

    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[3])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dxyc(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-omega covariance (dxy/c) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxyc values (one vector per jet)
"""
function get_dxyc(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[4])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_phic(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-omega covariance (phi/c) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phiomega values (one vector per jet)
"""
function get_phic(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[5])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dptdpt(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega covariance (dpt/dpt) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_omega_cov, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dptdpt values (one vector per jet)
"""
function get_dptdpt(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
                    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[6])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dxydz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-z0 covariance (dxy/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxy/dz values (one vector per jet)
"""
function get_dxydz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
                    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[7])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_phidz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-z0 covariance (dphi/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phidz values (one vector per jet)
"""
function get_phidz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[8])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_cdz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega-z0 covariance (c/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxdz values (one vector per jet)
"""
function get_cdz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[9])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dzdz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the z0 covariance (dz/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxdz values (one vector per jet)
"""
function get_dzdz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})

    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[10])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dxyctgtheta(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-tanLambda covariance (dxy/ctgtheta) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxyctgtheta values (one vector per jet)
"""
function get_dxyctgtheta(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[11])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_phictgtheta(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-tanLambda covariance (phi/ctgtheta) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phictgtheta values (one vector per jet)
"""
function get_phictgtheta(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[12])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_cctgtheta(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega-tanLambda covariance (c/ctgtheta) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of cctgtheta values (one vector per jet)
"""
function get_cctgtheta(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[13])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_dlambdadz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the tanLambda-z0 covariance (dlambda/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dlambdadz values (one vector per jet)
"""
function get_dlambdadz(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})

    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[14])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

"""
    get_detadeta(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the tanLambda covariance (deta/deta) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_tanlambda_cov, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of detadeta values (one vector per jet)
"""
function get_detadeta(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                    tracks::StructVector{EDM4hep.TrackState})

    result = Vector{JetConstituentsData}()
    tracks_len = length(tracks)

    for jc in jcs
        jet_result = Float32[]
        for p in jc
            if p.tracks.first < tracks_len
                push!(jet_result, tracks[p.tracks.first + 1].covMatrix[15])
            else
                push!(jet_result, -9.0f0)
            end
        end
        push!(result, jet_result)
    end
    return result
end

#################################################################################
# Clustered jet specific variables
#################################################################################


"""
    get_erel_log_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}) -> Vector{JetConstituentsData}

Calculate log of relative energy (log(E_const/E_jet)) for each constituent particle in clustered jets.
"""
function get_erel_log_cluster(jets::Vector{EEJet}, 
                            jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}})
    # Define the result type
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        jet_csts = Float32[]
        
        # Get jet energy
        e_jet = jets[i].E  # Assuming EEJet has an e property for energy
        
        # Get constituents for this jet
        if i <= length(jcs)
            constituents = jcs[i]
            
            for jc in constituents
                # Calculate relative energy and log
                val = (e_jet > 0.0) ? jc.energy / e_jet : 1.0
                erel_log = log10(val)
                push!(jet_csts, erel_log)
            end
        end
        
        push!(result, jet_csts)
    end
    
    return result
end

"""
    get_thetarel_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}) -> Vector{JetConstituentsData}

Calculate relative theta angle between constituent particle and clustered jet axis.
"""
function get_thetarel_cluster(jets::Vector{EEJet}, 
                            jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}})
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        jet_csts = Float32[]
        
        # Create a 4-vector for the jet
        jet = jets[i]
        px, py, pz, E = jet.px, jet.py, jet.pz, jet.E
        
        # Calculate jet direction angles
        p_mag = sqrt(px^2 + py^2 + pz^2)
        theta_jet = p_mag > 0 ? acos(pz / p_mag) : 0.0
        phi_jet = atan(py, px)
        
        # Get constituents for this jet
        if i <= length(jcs)
            constituents = jcs[i]
            
            for constituent in constituents
                # Create a 4-vector for the constituent
                p_const_x = constituent.momentum.x
                p_const_y = constituent.momentum.y
                p_const_z = constituent.momentum.z
                
                # Rotate the constituent vector to align with jet axis
                
                # First rotate around z-axis by -phi_jet
                p_rot_x = p_const_x * cos(-phi_jet) - p_const_y * sin(-phi_jet)
                p_rot_y = p_const_x * sin(-phi_jet) + p_const_y * cos(-phi_jet)
                p_rot_z = p_const_z
                
                # Then rotate around y-axis by -theta_jet
                p_rot2_x = p_rot_x * cos(-theta_jet) - p_rot_z * sin(-theta_jet)
                p_rot2_z = p_rot_x * sin(-theta_jet) + p_rot_z * cos(-theta_jet)
                p_rot2_y = p_rot_y
                
                # Calculate theta in rotated frame
                p_rot_mag = sqrt(p_rot2_x^2 + p_rot2_y^2 + p_rot2_z^2)
                theta_rel = p_rot_mag > 0 ? acos(p_rot2_z / p_rot_mag) : 0.0
                
                push!(jet_csts, theta_rel)
            end
        end
        
        push!(result, jet_csts)
    end
    
    return result
end

"""
    get_phirel_cluster(jets::Vector{EEJet}, 
                        jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}) -> Vector{JetConstituentsData}

Calculate relative phi angle between constituent particle and clustered jet axis.
"""
function get_phirel_cluster(jets::Vector{EEJet}, 
                            jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}})
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        jet_csts = Float32[]
        
        # Create a 4-vector for the jet
        jet = jets[i]
        # Access momentum components from EEJet
        px, py, pz, E = jet.px, jet.py, jet.pz, jet.E
        
        # Calculate jet direction angles
        p_mag = sqrt(px^2 + py^2 + pz^2)
        theta_jet = p_mag > 0 ? acos(pz / p_mag) : 0.0
        phi_jet = atan(py, px)
        
        # Get constituents for this jet
        if i <= length(jcs)
            constituents = jcs[i]
            
            for constituent in constituents
                # Get constituent momentum
                p_const_x = constituent.momentum.x
                p_const_y = constituent.momentum.y
                p_const_z = constituent.momentum.z
                
                # Rotate the constituent vector to align with jet axis
                
                # First rotate around z-axis by -phi_jet
                p_rot_x = p_const_x * cos(-phi_jet) - p_const_y * sin(-phi_jet)
                p_rot_y = p_const_x * sin(-phi_jet) + p_const_y * cos(-phi_jet)
                p_rot_z = p_const_z
                
                # Then rotate around y-axis by -theta_jet
                p_rot2_x = p_rot_x * cos(-theta_jet) - p_rot_z * sin(-theta_jet)
                p_rot2_z = p_rot_x * sin(-theta_jet) + p_rot_z * cos(-theta_jet)
                p_rot2_y = p_rot_y
                
                # Calculate phi in rotated frame
                phi_rel = atan(p_rot2_y, p_rot2_x)
                
                push!(jet_csts, phi_rel)
            end
        end
        
        push!(result, jet_csts)
    end
    
    return result
end

#################################################################################
# Particle ID variables
#################################################################################

"""
    get_isMu(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a muon.
"""
function get_isMu(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for jet_constituents in jcs
        is_mu = JetConstituentsData()
        
        for p in jet_constituents
            if abs(p.charge) > 0 && abs(p.mass - 0.105658) < 1e-3
                push!(is_mu, 1.0f0)
            else
                push!(is_mu, 0.0f0)
            end
        end
        
        push!(result, is_mu)
    end
    
    return result
end

"""
    get_isEl(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is an electron.
"""
function get_isEl(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for jet_constituents in jcs
        is_el = JetConstituentsData()
        
        for p in jet_constituents
            if abs(p.charge) > 0 && abs(p.mass - 0.000510999) < 1e-5
                push!(is_el, 1.0f0)
            else
                push!(is_el, 0.0f0)
            end
        end
        
        push!(result, is_el)
    end
    
    return result
end

"""
    get_isChargedHad(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a charged hadron.
"""
function get_isChargedHad(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for jet_constituents in jcs
        is_charged_had = JetConstituentsData()
        
        for p in jet_constituents
            if abs(p.charge) > 0 && abs(p.mass - 0.13957) < 1e-3
                push!(is_charged_had, 1.0f0)
            else
                push!(is_charged_had, 0.0f0)
            end
        end
        
        push!(result, is_charged_had)
    end
    
    return result
end

"""
    get_isGamma(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a photon (gamma).
"""
function get_isGamma(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for jet_constituents in jcs
        is_gamma = JetConstituentsData()
        
        for p in jet_constituents
            if p.type == 22  # PDG code for photon
                push!(is_gamma, 1.0f0)
            else
                push!(is_gamma, 0.0f0)
            end
        end
        
        push!(result, is_gamma)
    end
    
    return result
end

"""
    get_isNeutralHad(jcs::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a neutral hadron.
"""
function get_isNeutralHad(jcs::Vector{JetConstituents})
    result = Vector{JetConstituentsData}()
    
    for jet_constituents in jcs
        is_neutral_had = JetConstituentsData()
        
        for p in jet_constituents
            if p.type == 130  # PDG code for K_L^0 (common neutral hadron)
                push!(is_neutral_had, 1.0f0)
            else
                push!(is_neutral_had, 0.0f0)
            end
        end
        
        push!(result, is_neutral_had)
    end
    
    return result
end

#################################################################################
# Time-of-flight mass calculation
#################################################################################

"""
    get_mtof(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            track_L::AbstractArray{T} where T <: AbstractFloat,
            trackdata::StructVector{EDM4hep.Track},
            trackerhits::StructVector{EDM4hep.TrackerHit},
            gammadata::StructVector{EDM4hep.Cluster},
            nhdata::StructVector{EDM4hep.Cluster},
            calohits::StructVector{EDM4hep.CalorimeterHit},
            V::LorentzVector) -> Vector{JetConstituentsData}

Calculate the mass using time-of-flight measurements for each particle in each jet.
This is a Julia implementation of the C++ function get_mtof.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- track_L: Vector of track lengths
- trackdata: StructVector of Track objects
- trackerhits: StructVector of TrackerHit objects 
- gammadata: StructVector of photon Cluster objects
- nhdata: StructVector of neutral hadron Cluster objects
- calohits: StructVector of CalorimeterHit objects
- V: LorentzVector representing the primary vertex position and time

# Returns
Vector of vectors of mtof values (one vector per jet)
"""
function get_mtof(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                track_L::AbstractArray{T} where T <: AbstractFloat,
                trackdata::StructVector{EDM4hep.Track},
                trackerhits::StructVector{EDM4hep.TrackerHit},
                gammadata::StructVector{EDM4hep.Cluster},
                nhdata::StructVector{EDM4hep.Cluster},
                calohits::StructVector{EDM4hep.CalorimeterHit},
                V::LorentzVector)
    
    c_light = 2.99792458e8
    result = Vector{JetConstituentsData}()
    tracks_len = length(trackdata)
    gamma_len = length(gammadata)
    nh_len = length(nhdata)

    for i in eachindex(jcs)
        jc = jcs[i]
        tmp = JetConstituentsData()
        
        for j in eachindex(jc)
            # particle = jc[j]
            if jc[j].clusters.first < nh_len + gamma_len
                if jc[j].type == 130
                    hit_idx = (nhdata[jc[j].clusters.first+1 - gamma_len].hits.first+1)

                    T = calohits[hit_idx].time
                    X = calohits[hit_idx].position.x
                    Y = calohits[hit_idx].position.y
                    Z = calohits[hit_idx].position.z

                    tof = T 
                    L = sqrt((X - V.x)^2 + (Y - V.y)^2 + (Z - V.z)^2) * 0.001
                    beta = L/(tof * c_light)
                    E = jc[j].energy
                    if beta < 1.0 && beta > 0.0
                        push!(tmp, E * sqrt(1.0 - beta^2))
                    else
                        push!(tmp, 9.0f0)  # Invalid measurement
                    end
                elseif jc[j].type == 22
                    push!(tmp, 0.0f0)  # Photons have zero mass
                end
            end

            if jc[j].tracks.first < tracks_len
                if abs(jc[j].charge) > 0 && abs(jc[j].mass - 0.000510999) < 1e-5
                    push!(tmp, 0.000510999f0)  # Electron mass
                elseif abs(jc[j].charge) > 0 && abs(jc[j].mass - 0.105658) < 1e-3
                    push!(tmp, 0.105658f0)  # Muon mass
                else
                    Tin = V.t * 1e-3 / c_light 
                    Tout = trackerhits[trackdata[jc[j].tracks.first+1].trackerHits.last].time
                    tof = Tout - Tin

                    L = track_L[jc[j].tracks.first+1] * 0.001 
                    beta =  L/(tof * c_light)
                    p = sqrt(jc[j].momentum.x^2 + jc[j].momentum.y^2 + jc[j].momentum.z^2)
                    if beta < 1.0 && beta > 0.0
                        push!(tmp, p * sqrt(1.0 / (beta^2) - 1.0))
                    else
                        push!(tmp, 0.13957039f0)  # Default to pion mass
                    end
                end
            end
        end
        push!(result, tmp)
    end
    
    return result
end

##################################################################################
# dN/dx
##################################################################################

"""
    get_dndx(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
            dNdx::StructVector{EDM4hep.Quantity},
            trackdata::StructVector{EDM4hep.Track}, 
            JetsConstituents_isChargedHad::Vector{JetConstituentsData}) -> Vector{JetConstituentsData}

Calculate dE/dx or dN/dx for each charged hadron in jets. Neutrals, muons, and electrons are set to 0.
Only charged hadrons have valid dN/dx values.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- dNdx: StructVector of Quantity objects containing dN/dx measurements (EFlowTrack_2)
- trackdata: StructVector of Track objects (EFlowTrack)
- JetsConstituents_isChargedHad: Vector of vectors indicating which particles are charged hadrons

# Returns
Vector of vectors of dN/dx values (one vector per jet, one value per constituent)
"""
function get_dndx(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                dNdx::StructVector{EDM4hep.Quantity},
                trackdata::StructVector{EDM4hep.Track}, 
                JetsConstituents_isChargedHad::Vector{JetConstituentsData})
    
    result = Vector{JetConstituentsData}()
    tracks_len = length(trackdata)

    for i in 1:length(jcs)
        isChargedHad = JetsConstituents_isChargedHad[i]
        tmp = JetConstituentsData()
        for j in 1:length(jcs[i])
            if jcs[i][j].tracks.first + 1 <= tracks_len && Int(isChargedHad[j]) == 1
                push!(tmp, -1.0f0)
            else
                push!(tmp, 0.0f0)  # Not a charged hadron or no valid track
            end
        end
        push!(result, tmp)
    end

    return result
end

#################################################################################
# B-tagging variables, signed impact parameter, dist, etc.
#################################################################################

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

    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        p = Vector2(jets[i].px, jets[i].py)
        
        sip2d_values = JetConstituentsData()
        
        for j in eachindex(D0[i])
            if D0[i][j] != -9.0f0
                # Calculate the 2D impact point vector
                d0 = Vector2(-D0[i][j] * sin(phi0[i][j]), 
                              D0[i][j] * cos(phi0[i][j]))
                dot_product = d0.X * p.X + d0.Y * p.Y

                # Sign based on dot product and magnitude based on D0
                signed_ip = sign(dot_product) * abs(D0[i][j])
                
                push!(sip2d_values, signed_ip)
            else
                push!(sip2d_values, -9.0f0)
            end
        end
        
        push!(result, sip2d_values)
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
    
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(Sip2dVals)
        sig_values = JetConstituentsData()
        
        for j in eachindex(Sip2dVals[i])
            # Only calculate significance if the error is positive
            if j <= length(err2_D0[i]) && err2_D0[i][j] > 0.0
                # Calculate significance by dividing the value by its error
                significance = Sip2dVals[i][j] / sqrt(err2_D0[i][j])
                push!(sig_values, significance)
            else
                # Invalid measurement
                push!(sig_values, -9.0f0)
            end
        end
        
        push!(result, sig_values)
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
    result = Vector{JetConstituentsData}()
    for i in eachindex(jets)
        p = Vector3(jets[i].px, jets[i].py, jets[i].pz)
        cprojs = JetConstituentsData()
        
        for j in eachindex(D0[i])
            if D0[i][j] != -9.0
                # Create 3D vector of displacement at point of closest approach
                d = Vector3(-D0[i][j] * sin(phi0[i][j]), 
                            D0[i][j] * cos(phi0[i][j]), 
                            Z0[i][j])
                
                # Sign the impact parameter based on the dot product with jet direction
                # dot product in 3D
                dot_prod = d.X * p.X + d.Y * p.Y + d.Z * p.Z
                sign_val = dot_prod > 0.0 ? 1.0 : -1.0
                impact_val = sqrt(D0[i][j]^2 + Z0[i][j]^2)
                cprojs_val = sign_val * impact_val
                push!(cprojs, cprojs_val)
            else
                push!(cprojs, -9.0f0)
            end
        end
        
        push!(result, cprojs)
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
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(Sip3dVals)
        sigs = JetConstituentsData()  # Changed to initialize sigs as a JetConstituentsData
        
        for j in eachindex(Sip3dVals[i])
            if err2_D0[i][j] > 0.0
                sig = Sip3dVals[i][j] / sqrt(err2_D0[i][j] + err2_Z0[i][j])
                push!(sigs, sig)
            else
                push!(sigs, -9.0f0)
            end
        end
        
        push!(result, sigs)
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
                            jcs::Vector{JetConstituents},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the jet distance value for each particle, measuring the distance between
the point of closest approach and the jet axis.
"""
function get_JetDistVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                                jcs::Vector{JetConstituents},
                                D0::Vector{JetConstituentsData},
                                Z0::Vector{JetConstituentsData},
                                phi0::Vector{JetConstituentsData},
                                Bz::Float32)
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(jets)
        # p_jet = Vector3(jets[i].px, jets[i].py, jets[i].pz)
        p_jet = Vector3(jets[i].px, jets[i].py, jets[i].pz)
        tmp = JetConstituentsData()
        
        for j in eachindex(D0[i])
            if D0[i][j] != -9.0 && j <= length(jcs[i])
                # Create 3D vector of displacement at point of closest approach
                d = Vector3(-D0[i][j] * sin(phi0[i][j]), 
                            D0[i][j] * cos(phi0[i][j]), 
                            Z0[i][j])
                
                # Calculate particle momentum
                p_ct = Vector3(jcs[i][j].momentum.x, 
                                jcs[i][j].momentum.y, 
                                jcs[i][j].momentum.z)
                
                # Jet origin
                r_jet = Vector3(0.0, 0.0, 0.0)
                
                # Normal vector to plane containing particle and jet momenta
                n = v3_unit(v3_cross(p_ct, p_jet))

                
                # Distance is projection of displacement onto normal vector
                dist = n.X * (d.X - r_jet.X) + 
                       n.Y * (d.Y - r_jet.Y) + 
                       n.Z * (d.Z - r_jet.Z)
                
                push!(tmp, dist)
            else
                push!(tmp, -9.0f0)
            end
        end
        
        push!(result, tmp)
    end
    
    return result
end

function get_btagJetDistVal(jets::Vector{JetReconstruction.EEJet},
                            jcs::Vector{JetConstituents},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32)
    # Simply call the implementation functio
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
    result = Vector{JetConstituentsData}()
    
    for i in eachindex(JetDistVal)
        tmp = JetConstituentsData()
        
        for j in eachindex(JetDistVal[i])
            if err2_D0[i][j] > 0.0
                # 3D error
                err3d = sqrt(err2_D0[i][j] + err2_Z0[i][j])
                # Calculate significance
                jetdistsig = JetDistVal[i][j] / err3d
                push!(tmp, jetdistsig)
            else
                push!(tmp, -9.0f0)
            end
        end
        
        push!(result, tmp)
    end
    
    return result
end

function get_btagJetDistSig(JetDistVal::Vector{JetConstituentsData},
                            err2_D0::Vector{JetConstituentsData},
                            err2_Z0::Vector{JetConstituentsData})
    # Simply call the implementation function
    return get_JetDistSig(JetDistVal, err2_D0, err2_Z0)
end

#################################################################################
# Counting variables
#################################################################################

"""
    count_jets(jets::Vector{JetConstituents}) -> Int

Count the number of jets.
"""
function count_jets(jets::Vector{JetConstituents})
    return length(jets)
end

"""
    count_consts(jets::Vector{JetConstituents}) -> Vector{Int}

Count the number of constituents in each jet.
"""
function count_consts(jets::Vector{JetConstituents})
    result = Vector{Int}()
    
    for i in eachindex(jets)
        push!(result, length(jets[i]))
    end
    
    return result
end

"""
    count_type(isType::Vector{JetConstituentsData}) -> Vector{Int}

Count the number of particles of a specific type in each jet.
"""
function count_type(isType::Vector{JetConstituentsData})
    result = Vector{Int}()
    
    for i in eachindex(isType)
        count = 0
        for j in eachindex(isType[i])
            if Int(isType[i][j]) == 1
                count += 1
            end
        end
        push!(result, count)
    end
    
    return result
end
