using EDM4hep
using JetReconstruction
using StructArrays: StructVector

const JetConstituents = StructVector{ReconstructedParticle, <:Any}
const JetConstituentsData = Vector{Float32}

### Track Related Functions (5)

## Track Parameter Transformations (XPtoPar)
# XPtoPar_dxy - Transformed transverse impact parameter
# XPtoPar_dz - Transformed longitudinal impact parameter
# XPtoPar_phi - Transformed azimuthal angle
# XPtoPar_C - Track curvature parameter
# XPtoPar_ct - c×tau parameter

"""
    get_dxy(jcs::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the transverse impact parameter dxy for each particle in each jet relative to vertex V.
Reference: FCCAnalyses c++ function XPtoPar_dxy, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of dxy values (one vector per jet)
"""
function get_dxy(jcs::Vector{<:JetConstituents}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                V::LorentzVector, Bz::Float32)

    cSpeed_Bz = 2.99792458e8 * 1.0f-9 * Bz
    n_tracks = length(tracks)
    
    Vx, Vy = Float32(V.x), Float32(V.y)
    
    D0_values = tracks.D0
    phi_values = tracks.phi
    
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        charges = jet_constituents.charge
        track_indices = jet_constituents.tracks
        
        Float32[@inbounds(
            let track_idx = track_indices[i].first
                if track_idx < n_tracks
                    idx = track_idx + 1
                    D0 = D0_values[idx]
                    phi0 = phi_values[idx]
                    
                    sin_phi, cos_phi = sincos(phi0)
                    x1 = -D0 * sin_phi - Vx
                    x2 = D0 * cos_phi - Vy
                    
                    px = mom_x[i]
                    py = mom_y[i]
                    
                    a = -charges[i] * cSpeed_Bz
                    pt = sqrt(px^2 + py^2)
                    r2 = x1^2 + x2^2
                    cross = x1 * py - x2 * px
                    
                    # Compute impact parameter
                    discriminant = pt^2 - 2 * a * cross + a^2 * r2
                    if discriminant > 0
                        T = sqrt(discriminant)
                        if pt < 10.0f0
                            (T - pt) / a
                        else
                            (-2 * cross + a * r2) / (T + pt)
                        end
                    else
                        -9.0f0
                    end
                else
                    -9.0f0
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end


"""
    get_dz(jcs::Vector{JetConstituents}, 
           tracks::StructVector{EDM4hep.TrackState}, 
           V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the longitudinal impact parameter dz for each particle in each jet relative to vertex V.
Reference: FCCAnalyses c++ function XPtoPar_dz, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of dz values (one vector per jet)
"""
function get_dz(jcs::Vector{<:JetConstituents}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                V::LorentzVector, Bz::Float32)

    cSpeed_Bz = 2.99792458e8 * 1.0f-9 * Bz
    n_tracks = length(tracks)
    
    Vx, Vy, Vz = Float32(V.x), Float32(V.y), Float32(V.z)
    
    D0_values = tracks.D0
    Z0_values = tracks.Z0
    phi_values = tracks.phi
    
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        mom_z = jet_constituents.momentum.z
        charges = jet_constituents.charge
        track_indices = jet_constituents.tracks
        
        Float32[@inbounds(
            let track_idx = track_indices[i].first
                if track_idx < n_tracks
                    idx = track_idx + 1
                    D0 = D0_values[idx]
                    Z0 = Z0_values[idx]
                    phi0 = phi_values[idx]
                    
                    sin_phi, cos_phi = sincos(phi0)
                    
                    x1 = -D0 * sin_phi - Vx
                    x2 = D0 * cos_phi - Vy
                    x3 = Z0 - Vz
                    
                    px = mom_x[i]
                    py = mom_y[i]
                    pz = mom_z[i]
                    
                    # Compute intermediate values
                    a = -charges[i] * cSpeed_Bz
                    pt = sqrt(px^2 + py^2)
                    C = a / (2 * pt)
                    r2 = x1^2 + x2^2
                    cross = x1 * py - x2 * px
                    T = sqrt(pt^2 - 2 * a * cross + a^2 * r2)
                    
                    D = if pt < 10.0f0
                        (T - pt) / a
                    else
                        (-2 * cross + a * r2) / (T + pt)
                    end
                    
                    B_arg = max(r2 - D^2, 0.0f0) / (1 + 2 * C * D)
                    B = C * sqrt(B_arg)
                    if abs(B) > 1.0f0
                        B = sign(B)
                    end
                    
                    # Calculate st and ct
                    st = asin(B) / C
                    ct = pz / pt
                    
                    # Calculate z0
                    dot = x1 * px + x2 * py
                    if dot > 0.0f0
                        x3 - ct * st
                    else
                        x3 + ct * st
                    end
                else
                    -9.0f0
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_phi0(jcs::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the phi angle at the point of closest approach for each particle relative to vertex V.
This is a Julia implementation of the C++ function XPtoPar_phi.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of phi values (one vector per jet)
"""
function get_phi0(jcs::Vector{<:JetConstituents}, 
                  tracks::StructVector{EDM4hep.TrackState}, 
                  V::LorentzVector, Bz::Float32)

    cSpeed_Bz = 2.99792458e8 * 1.0f-9 * Bz
    n_tracks = length(tracks)
    
    Vx, Vy = Float32(V.x), Float32(V.y)
    
    D0_values = tracks.D0
    phi_values = tracks.phi
    
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        charges = jet_constituents.charge
        track_indices = jet_constituents.tracks
        
        Float32[@inbounds(
            let track_idx = track_indices[i].first
                if track_idx < n_tracks
                    idx = track_idx + 1
                    D0 = D0_values[idx]
                    phi0_track = phi_values[idx]

                    sin_phi, cos_phi = sincos(phi0_track)

                    x1 = -D0 * sin_phi - Vx
                    x2 = D0 * cos_phi - Vy
                    
                    px = mom_x[i]
                    py = mom_y[i]
                    
                    a = -charges[i] * cSpeed_Bz
                    
                    # Minimize redundant calculations
                    pt2 = px^2 + py^2
                    r2 = x1^2 + x2^2
                    cross = x1 * py - x2 * px
                    two_a_cross = 2 * a * cross
                    a2_r2 = a^2 * r2
                    
                    T = sqrt(pt2 - two_a_cross + a2_r2)
                    inv_T = 1.0f0 / T
                    
                    a_x1 = a * x1
                    a_x2 = a * x2
                    atan((py - a_x1) * inv_T, (px + a_x2) * inv_T)
                else
                    -9.0f0
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_c(jcs::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the track curvature for each particle in each jet.
Reference: FCCAnalyses c++ function XPtoPar_C, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of C values (one vector per jet)
"""
function get_c(jcs::Vector{<:JetConstituents}, 
                        tracks::StructVector{EDM4hep.TrackState}, 
                        Bz::Float32)

    cSpeed_Bz_half = 2.99792458e8 * 1.0f3 * 1.0f-15 * Bz * 0.5f0
    n_tracks = length(tracks)
    
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        charges = jet_constituents.charge
        track_indices = jet_constituents.tracks
        
        Float32[@inbounds(
            let track_idx = track_indices[i].first
                if track_idx < n_tracks
                    px = mom_x[i]
                    py = mom_y[i]
                    inv_pt = 1.0f0 / sqrt(px^2 + py^2)
                    copysign(cSpeed_Bz_half * inv_pt, charges[i])
                else
                    -9.0f0
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end

"""
    get_ct(jcs::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the c*tau for each particle in each jet.
Reference: FCCAnalyses c++ function XPtoPar_ct, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of ct values (one vector per jet)
"""
function get_ct(jcs::Vector{<:JetConstituents}, 
                tracks::StructVector{EDM4hep.TrackState}, 
                Bz::Float32)

    n_tracks = length(tracks)
    
    return [begin
        mom_x = jet_constituents.momentum.x
        mom_y = jet_constituents.momentum.y
        mom_z = jet_constituents.momentum.z
        track_indices = jet_constituents.tracks
        
        Float32[@inbounds(
            let track_idx = track_indices[i].first
                if track_idx < n_tracks
                    px = mom_x[i]
                    py = mom_y[i]
                    pz = mom_z[i]
                    pt = sqrt(px^2 + py^2)
                    pz / pt
                else
                    -9.0f0
                end
            end
        ) for i in eachindex(mom_x)]
    end for jet_constituents in jcs]
end