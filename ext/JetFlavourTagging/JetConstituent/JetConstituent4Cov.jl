using EDM4hep
using JetReconstruction
using StructArrays: StructVector

### Covariance Matrix Elements (15)

## Diagonal Elements
# get_omega_cov - Omega variance
# get_d0_cov - d0 variance
# get_z0_cov - z0 variance
# get_phi0_cov - phi0 variance
# get_tanlambda_cov - tanLambda variance

## Off-diagonal Elements
# get_d0_z0_cov
# get_phi0_d0_cov
# get_phi0_z0_cov
# get_tanlambda_phi0_cov
# get_tanlambda_d0_cov
# get_tanlambda_z0_cov
# get_omega_tanlambda_cov
# get_omega_phi0_cov
# get_omega_d0_cov
# get_omega_z0_cov

"""
    get_dxydxy(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0 covariance (dxy/dxy) for each particle.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxydxy values
"""
function get_dxydxy(jcs::Vector{<:JetConstituents}, 
                    tracks::StructVector{EDM4hep.TrackState})
    
    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[1], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dphidxy(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-d0 covariance (dphi/dxy) for each particle.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dphidxy values 
"""
function get_dphidxy(jcs::Vector{<:JetConstituents}, 
                     tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[2], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dphidphi(jcs::Vector{JetConstituents}, 
                tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi covariance (dphi/dphi) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_phi0_cov, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dphidphi values 
"""
function get_dphidphi(jcs::Vector{<:JetConstituents}, 
                      tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[3], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dxyc(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-omega covariance (dxy/c) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxyc values 
"""
function get_dxyc(jcs::Vector{<:JetConstituents}, 
                  tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[4], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_phic(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-omega covariance (phi/c) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phiomega values
"""
function get_phic(jcs::Vector{<:JetConstituents}, 
                  tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[5], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dptdpt(jcs::Vector{JetConstituents}, 
               tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega covariance (dpt/dpt) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_omega_cov, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dptdpt values
"""
function get_dptdpt(jcs::Vector{<:JetConstituents}, 
                    tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[6], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dxydz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-z0 covariance (dxy/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxy/dz values
"""
function get_dxydz(jcs::Vector{<:JetConstituents}, 
                   tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[7], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_phidz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-z0 covariance (dphi/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phidz values
"""
function get_phidz(jcs::Vector{<:JetConstituents}, 
                   tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[8], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_cdz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega-z0 covariance (c/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxdz values
"""
function get_cdz(jcs::Vector{<:JetConstituents}, 
                 tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[9], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dzdz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the z0 covariance (dz/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of z0 covariance values
"""
function get_dzdz(jcs::Vector{<:JetConstituents}, 
                  tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[10], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dxyctgtheta(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-tanLambda covariance (dxy/ctgtheta) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxyctgtheta values
"""
function get_dxyctgtheta(jcs::Vector{<:JetConstituents}, 
                         tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[11], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_phictgtheta(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-tanLambda covariance (phi/ctgtheta) for each particle.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phictgtheta values 
"""
function get_phictgtheta(jcs::Vector{JetConstituents}, 
                         tracks::StructVector{EDM4hep.TrackState})
    
    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[12], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_cctgtheta(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega-tanLambda covariance (c/ctgtheta) for each particle.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of cctgtheta values
"""
function get_cctgtheta(jcs::Vector{JetConstituents}, 
                       tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[13], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_dlambdadz(jcs::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the tanLambda-z0 covariance (dlambda/dz) for each particle.

# Arguments
- jcs: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dlambdadz values 
"""
function get_dlambdadz(jcs::Vector{JetConstituents}, 
                       tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[14], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end

"""
    get_detadeta(jcs::Vector{JetConstituents}, 
                 tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the tanLambda covariance (deta/deta) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_tanlambda_cov, adapted for jet constituents.

# Arguments
- jcs: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of detadeta values
"""
function get_detadeta(jcs::Vector{JetConstituents}, 
                      tracks::StructVector{EDM4hep.TrackState})

    n_jets = length(jcs)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)
    
    @inbounds for j in 1:n_jets
        jet_constituents = jcs[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)
        
        jet_result = Vector{Float32}(undef, n_particles)
        
        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks, 
                                   tracks[track_idx + 1].covMatrix[15], 
                                   -9.0f0)
        end
        
        result[j] = jet_result
    end
    
    return result
end