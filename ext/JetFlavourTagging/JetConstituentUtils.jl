using EDM4hep
using JetReconstruction
using StructArrays: StructVector

# Import physical constants
include("JetPhysicalConstants.jl")
using .JetPhysicalConstants

const JetConstituents = StructVector{ReconstructedParticle, <:Any}
const JetConstituentsData = Vector{Float32}

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

"""
    get_pt(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the transverse momentum of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents (each element contains particles for one jet)

# Returns
A vector of vectors of transverse momentum values (sqrt(px^2 + py^2))
"""
function get_pt(jets_constituents::Vector{<:JetConstituents})
    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                Float32[@inbounds sqrt(mom_x[i]^2 + mom_y[i]^2)
                        for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_p(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the momentum magnitude of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of momentum magnitudes (sqrt(px^2 + py^2 + pz^2))
"""
function get_p(jets_constituents::Vector{<:JetConstituents})
    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                mom_z = jet_constituents.momentum.z
                Float32[@inbounds sqrt(mom_x[i]^2 + mom_y[i]^2 + mom_z[i]^2)
                        for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_e(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the energy of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of energy values
"""
function get_e(jets_constituents::Vector{<:JetConstituents})
    return [jet_constituents.energy for jet_constituents in jets_constituents]
end

"""
    get_type(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the PDG type of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of particle types (PDG codes/Particle IDs)
"""
function get_type(jets_constituents::Vector{<:JetConstituents})
    return [jet_constituents.type for jet_constituents in jets_constituents]
end

"""
    get_mass(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the mass of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of mass values
"""
function get_mass(jets_constituents::Vector{<:JetConstituents})
    return [jet_constituents.mass for jet_constituents in jets_constituents]
end

"""
    get_charge(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the charge of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of charge values
"""
function get_charge(jets_constituents::Vector{<:JetConstituents})
    return [jet_constituents.charge for jet_constituents in jets_constituents]
end

"""
    get_theta(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the polar angle of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of polar angle values
"""
function get_theta(jets_constituents::Vector{<:JetConstituents})
    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                mom_z = jet_constituents.momentum.z
                Float32[@inbounds(let x = mom_x[i], y = mom_y[i], z = mom_z[i]
                                      (x == 0.0f0 && y == 0.0f0 && z == 0.0f0) ? 0.0f0 :
                                      atan(sqrt(x^2 + y^2), z)
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_phi(jets_constituents::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Get the azimuthal angle of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of azimuthal angle values
"""
function get_phi(jets_constituents::Vector{<:JetConstituents})
    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                Float32[@inbounds(let x = mom_x[i], y = mom_y[i]
                                      (x == 0.0f0 && y == 0.0f0) ? 0.0f0 : atan(y, x)
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_y(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}
    
Get the rapidity of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of rapidity values
"""
function get_y(jets_constituents::Vector{<:JetConstituents})
    return [begin
                energies = jet_constituents.energy
                mom_z = jet_constituents.momentum.z
                Float32[@inbounds(let e = energies[i], pz = mom_z[i]
                                      0.5f0 * log((e + pz) / (e - pz))
                                  end) for i in eachindex(energies)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_eta(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Get the pseudorapidity of each particle in each jet.

# Arguments
- jets_constituents: Vector of jet constituents

# Returns
A vector of vectors of pseudorapidity values (eta = -ln(tan(theta/2)))
"""
function get_eta(jets_constituents::Vector{<:JetConstituents})
    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                mom_z = jet_constituents.momentum.z
                Float32[@inbounds(let x = mom_x[i], y = mom_y[i], z = mom_z[i]
                                      p = sqrt(x^2 + y^2 + z^2)
                                      if p == 0.0f0
                                          0.0f0
                                      elseif p == abs(z)  # particle along beam axis
                                          sign(z) * Inf32
                                      else
                                          0.5f0 * log((p + z) / (p - z))
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_Bz(jets_constituents::Vector{<:JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Calculate the magnetic field Bz for each particle based on track curvature and momentum.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: Vector of track states (used to get the omega value)

# Returns 
A vector of vectors of Bz values.
"""
function get_Bz(jets_constituents::Vector{<:JetConstituents},
                tracks::StructVector{EDM4hep.TrackState})
    a = C_LIGHT * MM_TO_M / FS_TO_S
    n_tracks = length(tracks)

    # If tracks is a StructVector, we can access omega column directly
    omega_values = tracks.omega

    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                charges = jet_constituents.charge
                track_indices = jet_constituents.tracks

                Float32[@inbounds(let track_idx = track_indices[i].first
                                      if track_idx < n_tracks
                                          pt = sqrt(mom_x[i]^2 + mom_y[i]^2)
                                          omega_values[track_idx + 1] / a * pt *
                                          copysign(1.0f0, charges[i])
                                      else
                                          UNDEF_VAL
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

### Track Related Functions (5)

## Track Parameter Transformations (XPtoPar)
# XPtoPar_dxy - Transformed transverse impact parameter
# XPtoPar_dz - Transformed longitudinal impact parameter
# XPtoPar_phi - Transformed azimuthal angle
# XPtoPar_C - Track curvature parameter
# XPtoPar_ct - c×tau parameter

"""
    get_dxy(jets_constituents::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the transverse impact parameter dxy for each particle in each jet relative to vertex V.
Reference: FCCAnalyses c++ function XPtoPar_dxy, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of dxy values (one vector per jet)
"""
function get_dxy(jets_constituents::Vector{<:JetConstituents},
                 tracks::StructVector{EDM4hep.TrackState},
                 V::LorentzVector, Bz::Float32)
    cSpeed_Bz = C_LIGHT * NS_TO_S * Bz
    n_tracks = length(tracks)

    Vx, Vy = Float32(V.x), Float32(V.y)

    D0_values = tracks.D0
    phi_values = tracks.phi

    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                charges = jet_constituents.charge
                track_indices = jet_constituents.tracks

                Float32[@inbounds(let track_idx = track_indices[i].first
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
                                              t = sqrt(discriminant)
                                              if pt < 10.0f0
                                                  (t - pt) / a
                                              else
                                                  (-2 * cross + a * r2) / (t + pt)
                                              end
                                          else
                                              UNDEF_VAL
                                          end
                                      else
                                          UNDEF_VAL
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_dz(jets_constituents::Vector{JetConstituents}, 
           tracks::StructVector{EDM4hep.TrackState}, 
           V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the longitudinal impact parameter dz for each particle in each jet relative to vertex V.
Reference: FCCAnalyses c++ function XPtoPar_dz, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of dz values (one vector per jet)
"""
function get_dz(jets_constituents::Vector{<:JetConstituents},
                tracks::StructVector{EDM4hep.TrackState},
                V::LorentzVector, Bz::Float32)
    cSpeed_Bz = C_LIGHT * NS_TO_S * Bz
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

                Float32[@inbounds(let track_idx = track_indices[i].first
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
                                          c = a / (2 * pt)
                                          r2 = x1^2 + x2^2
                                          cross = x1 * py - x2 * px
                                          t = sqrt(pt^2 - 2 * a * cross + a^2 * r2)

                                          d = if pt < 10.0f0
                                              (t - pt) / a
                                          else
                                              (-2 * cross + a * r2) / (t + pt)
                                          end

                                          b_arg = max(r2 - d^2, 0.0f0) / (1 + 2 * c * d)
                                          b = c * sqrt(b_arg)
                                          if abs(b) > 1.0f0
                                              b = sign(b)
                                          end

                                          # Calculate st and ct
                                          st = asin(b) / c
                                          ct = pz / pt

                                          # Calculate z0
                                          dot = x1 * px + x2 * py
                                          if dot > 0.0f0
                                              x3 - ct * st
                                          else
                                              x3 + ct * st
                                          end
                                      else
                                          UNDEF_VAL
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_phi0(jets_constituents::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            V::LorentzVector, Bz::Float32) -> Vector{JetConstituentsData}

Calculate the phi angle at the point of closest approach for each particle relative to vertex V.
This is a Julia implementation of the C++ function XPtoPar_phi.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects
- V: LorentzVector representing the primary vertex
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of phi values (one vector per jet)
"""
function get_phi0(jets_constituents::Vector{<:JetConstituents},
                  tracks::StructVector{EDM4hep.TrackState},
                  V::LorentzVector, Bz::Float32)
    cSpeed_Bz = C_LIGHT * NS_TO_S * Bz
    n_tracks = length(tracks)

    Vx, Vy = Float32(V.x), Float32(V.y)

    D0_values = tracks.D0
    phi_values = tracks.phi

    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                charges = jet_constituents.charge
                track_indices = jet_constituents.tracks

                Float32[@inbounds(let track_idx = track_indices[i].first
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

                                          t = sqrt(pt2 - two_a_cross + a2_r2)
                                          inv_t = 1.0f0 / t

                                          a_x1 = a * x1
                                          a_x2 = a * x2
                                          atan((py - a_x1) * inv_t, (px + a_x2) * inv_t)
                                      else
                                          UNDEF_VAL
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_c(jets_constituents::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the track curvature for each particle in each jet.
Reference: FCCAnalyses c++ function XPtoPar_C, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of C values (one vector per jet)
"""
function get_c(jets_constituents::Vector{<:JetConstituents},
               tracks::StructVector{EDM4hep.TrackState},
               Bz::Float32)
    cSpeed_Bz_half = C_LIGHT * MM_TO_M / FS_TO_S * Bz * 0.5f0
    n_tracks = length(tracks)

    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                charges = jet_constituents.charge
                track_indices = jet_constituents.tracks

                Float32[@inbounds(let track_idx = track_indices[i].first
                                      if track_idx < n_tracks
                                          px = mom_x[i]
                                          py = mom_y[i]
                                          inv_pt = 1.0f0 / sqrt(px^2 + py^2)
                                          copysign(cSpeed_Bz_half * inv_pt, charges[i])
                                      else
                                          UNDEF_VAL
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

"""
    get_ct(jets_constituents::Vector{JetConstituents}, 
            tracks::StructVector{EDM4hep.TrackState}, 
            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the c*tau for each particle in each jet.
Reference: FCCAnalyses c++ function XPtoPar_ct, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects
- Bz: The magnetic field in Tesla

# Returns
Vector of vectors of ct values (one vector per jet)
"""
function get_ct(jets_constituents::Vector{<:JetConstituents},
                tracks::StructVector{EDM4hep.TrackState},
                Bz::Float32)
    n_tracks = length(tracks)

    return [begin
                mom_x = jet_constituents.momentum.x
                mom_y = jet_constituents.momentum.y
                mom_z = jet_constituents.momentum.z
                track_indices = jet_constituents.tracks

                Float32[@inbounds(let track_idx = track_indices[i].first
                                      if track_idx < n_tracks
                                          px = mom_x[i]
                                          py = mom_y[i]
                                          pz = mom_z[i]
                                          pt = sqrt(px^2 + py^2)
                                          pz / pt
                                      else
                                          UNDEF_VAL
                                      end
                                  end) for i in eachindex(mom_x)]
            end
            for jet_constituents in jets_constituents]
end

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
    get_dxydxy(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0 covariance (dxy/dxy) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxydxy values
"""
function get_dxydxy(jets_constituents::Vector{<:JetConstituents},
                    tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[1],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dphidxy(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-d0 covariance (dphi/dxy) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dphidxy values 
"""
function get_dphidxy(jets_constituents::Vector{<:JetConstituents},
                     tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[2],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dphidphi(jets_constituents::Vector{JetConstituents}, 
                tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi covariance (dphi/dphi) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_phi0_cov, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dphidphi values 
"""
function get_dphidphi(jets_constituents::Vector{<:JetConstituents},
                      tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[3],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dxyc(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-omega covariance (dxy/c) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxyc values 
"""
function get_dxyc(jets_constituents::Vector{<:JetConstituents},
                  tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[4],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_phic(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-omega covariance (phi/c) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phiomega values
"""
function get_phic(jets_constituents::Vector{<:JetConstituents},
                  tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[5],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dptdpt(jets_constituents::Vector{JetConstituents}, 
               tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega covariance (dpt/dpt) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_omega_cov, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents (each element contains particles for one jet)
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dptdpt values
"""
function get_dptdpt(jets_constituents::Vector{<:JetConstituents},
                    tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[6],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dxydz(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-z0 covariance (dxy/dz) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxy/dz values
"""
function get_dxydz(jets_constituents::Vector{<:JetConstituents},
                   tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[7],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_phidz(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-z0 covariance (dphi/dz) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phidz values
"""
function get_phidz(jets_constituents::Vector{<:JetConstituents},
                   tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[8],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_cdz(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega-z0 covariance (c/dz) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxdz values
"""
function get_cdz(jets_constituents::Vector{<:JetConstituents},
                 tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[9],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dzdz(jets_constituents::Vector{JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the z0 covariance (dz/dz) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of z0 covariance values
"""
function get_dzdz(jets_constituents::Vector{<:JetConstituents},
                  tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[10],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dxyctgtheta(jets_constituents::Vector{<:JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the d0-tanLambda covariance (dxy/ctgtheta) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dxyctgtheta values
"""
function get_dxyctgtheta(jets_constituents::Vector{<:JetConstituents},
                         tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[11],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_phictgtheta(jets_constituents::Vector{<:JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the phi0-tanLambda covariance (phi/ctgtheta) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of phictgtheta values 
"""
function get_phictgtheta(jets_constituents::Vector{<:JetConstituents},
                         tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[12],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_cctgtheta(jets_constituents::Vector{<:JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the omega-tanLambda covariance (c/ctgtheta) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of cctgtheta values
"""
function get_cctgtheta(jets_constituents::Vector{<:JetConstituents},
                       tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[13],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_dlambdadz(jets_constituents::Vector{<:JetConstituents}, tracks::Vector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the tanLambda-z0 covariance (dlambda/dz) for each particle.

# Arguments
- jets_constituents: Vector of jet constituents 
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of dlambdadz values 
"""
function get_dlambdadz(jets_constituents::Vector{<:JetConstituents},
                       tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[14],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

"""
    get_detadeta(jets_constituents::Vector{<:JetConstituents}, 
                 tracks::StructVector{EDM4hep.TrackState}) -> Vector{JetConstituentsData}

Get the tanLambda covariance (deta/deta) for each particle in each jet from its associated track.
Reference: FCCAnalyses c++ function get_tanlambda_cov, adapted for jet constituents.

# Arguments
- jets_constituents: Vector of jet constituents
- tracks: StructVector of TrackState objects

# Returns
Vector of vectors of detadeta values
"""
function get_detadeta(jets_constituents::Vector{<:JetConstituents},
                      tracks::StructVector{EDM4hep.TrackState})
    n_jets = length(jets_constituents)
    n_tracks = length(tracks)
    result = Vector{JetConstituentsData}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        track_indices = jet_constituents.tracks
        n_particles = length(track_indices)

        jet_result = Vector{Float32}(undef, n_particles)

        @simd ivdep for i in 1:n_particles
            track_idx = track_indices[i].first
            jet_result[i] = ifelse(track_idx < n_tracks,
                                   tracks[track_idx + 1].covMatrix[15],
                                   UNDEF_VAL)
        end

        result[j] = jet_result
    end

    return result
end

### Particle Identification (5)

# get_is_mu - Check if constituent is muon
# get_is_el - Check if constituent is electron
# get_is_charged_had - Check if constituent is charged hadron
# get_is_gamma - Check if constituent is photon
# get_is_neutral_had - Check if constituent is neutral hadron

"""
    get_is_mu(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a muon.

# Arguments
- jets_constituents: Vector of jet constituents 

# Returns
Vector of vectors of is muon boolean values as Float32.
"""
function get_is_mu(jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        charges = jet_constituents.charge
        masses = jet_constituents.mass
        n_particles = length(charges)

        is_mu = Vector{Float32}(undef, n_particles)

        @simd for i in 1:n_particles
            charge_check = abs(charges[i]) > 0
            mass_check = abs(masses[i] - MUON_MASS) < MUON_TOLERANCE
            is_mu[i] = (charge_check & mass_check) ? 1.0f0 : 0.0f0
        end

        result[j] = is_mu
    end

    return result
end

"""
    get_is_el(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is an electron.

# Arguments
- jets_constituents: Vector of jet constituents 

# Returns
Vector of vectors of is electron boolean values as Float32.
"""
function get_is_el(jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        charges = jet_constituents.charge
        masses = jet_constituents.mass
        n_particles = length(charges)

        is_mu = Vector{Float32}(undef, n_particles)

        @simd for i in 1:n_particles
            charge_check = abs(charges[i]) > 0
            mass_check = abs(masses[i] - ELECTRON_MASS) < ELECTRON_TOLERANCE
            is_mu[i] = (charge_check & mass_check) ? 1.0f0 : 0.0f0
        end

        result[j] = is_mu
    end

    return result
end

"""
    get_is_charged_had(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a charged hadron.
    
# Arguments
- jets_constituents: Vector of jet constituents 

# Returns
Vector of vectors of is charged hadron boolean values as Float32.
"""
function get_is_charged_had(jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        charges = jet_constituents.charge
        masses = jet_constituents.mass
        n_particles = length(charges)

        is_mu = Vector{Float32}(undef, n_particles)

        @simd for i in 1:n_particles
            charge_check = abs(charges[i]) > 0
            mass_check = abs(masses[i] - PION_MASS) < PION_TOLERANCE
            is_mu[i] = (charge_check & mass_check) ? 1.0f0 : 0.0f0
        end

        result[j] = is_mu
    end

    return result
end

"""
    get_is_gamma(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a photon (gamma) (PDG 22).

# Arguments
- jets_constituents: Vector of jet constituents 

# Returns
Vector of vectors of is photon boolean values as Float32.
"""
function get_is_gamma(jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        types = jet_constituents.type
        n_particles = length(types)

        is_gamma = Vector{Float32}(undef, n_particles)

        @simd for i in 1:n_particles
            is_gamma[i] = types[i] == PDG_PHOTON ? 1.0f0 : 0.0f0
        end

        result[j] = is_gamma
    end

    return result
end

"""
    get_is_neutral_had(jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Check if each constituent particle is a neutral hadron (PDG 130).

# Arguments
- jets_constituents: Vector of jet constituents 

# Returns
Vector of vectors of is neutral hadron boolean values as Float32.
"""
function get_is_neutral_had(jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for j in 1:n_jets
        jet_constituents = jets_constituents[j]
        types = jet_constituents.type
        n_particles = length(types)

        is_gamma = Vector{Float32}(undef, n_particles)

        @simd for i in 1:n_particles
            is_gamma[i] = types[i] == PDG_K_LONG ? 1.0f0 : 0.0f0
        end

        result[j] = is_gamma
    end

    return result
end

### Relative Kinematics (5)

# get_erel_cluster - Relative energy for clustered jets
# get_erel_log_cluster - Log of relative energy for clustered jets
# get_thetarel_cluster - Relative polar angle for clustered jets
# get_phirel_cluster - Relative azimuthal angle for clustered jets
# get_theta_phi_rel_cluster - Combined relative angles for clustered jets as they use the same logic

"""
    get_erel_cluster(jets::Vector{EEJet}, 
                        jets_constituents::Vector{<:JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative energy (E_const/E_jet) for each constituent particle in clustered jets.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets.
- `jets_constituents::Vector{<:JetConstituents}`: Vector of jet constituents corresponding to the jets.

# Returns
Vector containing relative energy values for each constituent in the jets.
"""
function get_erel_cluster(jets::Vector{EEJet},
                          jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets)
    res = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for i in 1:n_jets
        e_jet = jets[i].E
        constituents_collection = jets_constituents[i]
        energies = constituents_collection.energy
        n_constituents = length(energies)
        jet_constituents_collection = Vector{Float32}(undef, n_constituents)

        if e_jet > 0.0f0
            inv_e_jet = 1.0f0 / e_jet
            @inbounds @simd for j in 1:n_constituents
                jet_constituents_collection[j] = energies[j] * inv_e_jet
            end
        else
            @inbounds @simd for j in 1:n_constituents
                jet_constituents_collection[j] = 0.0f0
            end
        end

        res[i] = jet_constituents_collection
    end

    return res
end

"""
    get_erel_log_cluster(jets::Vector{EEJet}, 
                        jets_constituents::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate log of relative energy (log(E_const/E_jet)) for each constituent particle in clustered jets.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jets_constituents::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Vector containing log of relative energy values for each constituent in the jets
"""
function get_erel_log_cluster(jets::Vector{EEJet},
                              jets_constituents::Vector{<:JetConstituents})
    n_jets = length(jets)
    res = Vector{Vector{Float32}}(undef, n_jets)

    @inbounds for i in 1:n_jets
        e_jet = jets[i].E
        constituents_collection = jets_constituents[i]
        energies = constituents_collection.energy
        n_constituents = length(energies)
        jet_constituents_collection = Vector{Float32}(undef, n_constituents)

        if e_jet > 0.0f0
            inv_e_jet = 1.0f0 / e_jet
            @inbounds @simd for j in 1:n_constituents
                jet_constituents_collection[j] = log10(energies[j] * inv_e_jet)
            end
        else
            @inbounds @simd for j in 1:n_constituents
                jet_constituents_collection[j] = 0.0f0
            end
        end

        res[i] = jet_constituents_collection
    end

    return res
end

"""
    get_thetarel_cluster(jets::Vector{EEJet}, 
                        jets_constituents::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative theta angle between constituent particle and clustered jet axis.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jets_constituents::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Vector containing relative theta angle values for each constituent in the jets
"""
function get_thetarel_cluster(jets::Vector{EEJet},
                              jets_constituents::Vector{<:JetConstituents})
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

        constituents_collection = jets_constituents[i]
        mom_x = constituents_collection.momentum.x
        mom_y = constituents_collection.momentum.y
        mom_z = constituents_collection.momentum.z
        n_constituents = length(mom_x)
        jet_constituents_collection = Vector{Float32}(undef, n_constituents)

        @inbounds for j in 1:n_constituents
            # First rotation
            p_rot_x = mom_x[j] * cos_phi - mom_y[j] * sin_phi
            p_rot_y = mom_x[j] * sin_phi + mom_y[j] * cos_phi

            # Second rotation
            p_rot2_x = p_rot_x * cos_theta - mom_z[j] * sin_theta
            p_rot2_z = p_rot_x * sin_theta + mom_z[j] * cos_theta

            pt_rot_sq = p_rot2_x^2 + p_rot_y^2
            pt_rot = sqrt(pt_rot_sq)
            jet_constituents_collection[j] = atan(pt_rot, p_rot2_z)
        end

        result[i] = jet_constituents_collection
    end

    return result
end

"""
    get_phirel_cluster(jets::Vector{EEJet}, 
                        jets_constituents::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative phi angle between constituent particle and clustered jet axis.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jets_constituents::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Vector containing relative phi angle values for each constituent in the jets
"""
function get_phirel_cluster(jets::Vector{EEJet},
                            jets_constituents::Vector{<:JetConstituents})
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

        constituents_collection = jets_constituents[i]
        mom_x = constituents_collection.momentum.x
        mom_y = constituents_collection.momentum.y
        mom_z = constituents_collection.momentum.z
        n_constituents = length(mom_x)
        jet_constituents_collection = Vector{Float32}(undef, n_constituents)

        @inbounds for j in 1:n_constituents
            # First rotation around z-axis by -phi_jet
            p_rot_x = mom_x[j] * cos_phi - mom_y[j] * sin_phi
            p_rot_y = mom_x[j] * sin_phi + mom_y[j] * cos_phi

            # Second rotation around y-axis by -theta_jet
            p_rot2_x = p_rot_x * cos_theta - mom_z[j] * sin_theta

            # Calculate phi in rotated frame
            jet_constituents_collection[j] = atan(p_rot_y, p_rot2_x)
        end

        result[i] = jet_constituents_collection
    end

    return result
end

"""
    get_thetaphirel_cluster(jets::Vector{EEJet}, 
                        jets_constituents::Vector{JetConstituents}) -> Vector{JetConstituentsData}

Calculate relative theta and phi angles between constituent particles and clustered jet axis.

# Arguments
- `jets::Vector{EEJet}`: Vector of clustered jets
- `jets_constituents::Vector{JetConstituents}`: Vector of jet constituents corresponding to the jets

# Returns
Tuple of Vectors containing relative theta and phi angle values for each constituent in the jets
"""
function get_thetarel_phirel_cluster(jets::Vector{EEJet},
                                     jets_constituents::Vector{<:JetConstituents})
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

        constituents_collection = jets_constituents[i]
        mom_x = constituents_collection.momentum.x
        mom_y = constituents_collection.momentum.y
        mom_z = constituents_collection.momentum.z
        n_constituents = length(mom_x)

        jet_theta = Vector{Float32}(undef, n_constituents)
        jet_phi = Vector{Float32}(undef, n_constituents)

        @inbounds for j in 1:n_constituents
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

### Special Measurements (2)

# get_dndx - dE/dx measurement (energy loss)
# get_mtof - Mass from time-of-flight measurement

"""
    get_dndx(jets_constituents::Vector{JetConstituents}, 
            dNdx::StructVector{EDM4hep.Quantity},
            trackdata::StructVector{EDM4hep.Track}, 
            JetsConstituents_isChargedHad::Vector{JetConstituentsData}) -> Vector{JetConstituentsData}

Calculate dE/dx or dN/dx for each charged hadron in jets. Neutrals, muons, and electrons are set to 0.
Only charged hadrons have valid dN/dx values.

# Arguments
- jets_constituents: Vector of jet constituents (each element contains particles for one jet)
- dNdx: StructVector of Quantity objects containing dN/dx measurements (EFlowTrack_2)
- trackdata: StructVector of Track objects (EFlowTrack)
- JetsConstituents_isChargedHad: Vector of vectors indicating which particles are charged hadrons

# Returns
Vector of vectors of dN/dx values (one vector per jet, one value per constituent)
"""
function get_dndx(jets_constituents::Vector{<:JetConstituents},
                  dNdx::StructVector{EDM4hep.Quantity},
                  trackdata::StructVector{EDM4hep.Track},
                  JetsConstituents_isChargedHad::Vector{Vector{Float32}})
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)
    tracks_len = length(trackdata)

    @inbounds for i in 1:n_jets
        jet = jets_constituents[i]
        tracks_first = jet.tracks
        isChargedHad = JetsConstituents_isChargedHad[i]
        n_constituents = length(jet)
        tmp = Vector{Float32}(undef, n_constituents)

        @simd ivdep for j in 1:n_constituents
            has_valid_track = tracks_first[j].first + 1 <= tracks_len
            is_charged_had = isChargedHad[j] == 1.0f0
            tmp[j] = ifelse(has_valid_track & is_charged_had, -1.0f0, 0.0f0)
        end

        result[i] = tmp
    end

    return result
end

function get_mtof(jets_constituents::Vector{<:JetConstituents},
                  track_L::AbstractArray{T} where {T <: Float32},
                  trackdata::StructVector{EDM4hep.Track},
                  trackerhits::StructVector{EDM4hep.TrackerHit},
                  gammadata::StructVector{EDM4hep.Cluster},
                  nhdata::StructVector{EDM4hep.Cluster},
                  calohits::StructVector{EDM4hep.CalorimeterHit},
                  V::LorentzVector)
    n_jets = length(jets_constituents)
    result = Vector{Vector{Float32}}(undef, n_jets)

    # Pre-compute limits
    tracks_len = length(trackdata)
    gamma_len = length(gammadata)
    nh_len = length(nhdata)
    cluster_limit = nh_len + gamma_len

    # Pre-compute vertex values
    vx, vy, vz = V.x, V.y, V.z
    v_t_scaled = V.t * MM_TO_M * C_LIGHT_INV  # Tin calculation

    @inbounds for i in 1:n_jets
        single_jet_constituents = jets_constituents[i]
        n_constituents = length(single_jet_constituents)

        # Pre-allocate for this jet
        tmp = Vector{Float32}(undef, n_constituents)
        result_idx = 0

        # Access fields once
        clusters_first = single_jet_constituents.clusters
        tracks_first = single_jet_constituents.tracks
        types = single_jet_constituents.type
        charges = single_jet_constituents.charge
        masses = single_jet_constituents.mass
        energies = single_jet_constituents.energy
        mom_x = single_jet_constituents.momentum.x
        mom_y = single_jet_constituents.momentum.y
        mom_z = single_jet_constituents.momentum.z

        for j in 1:n_constituents
            cluster_idx = clusters_first[j].first
            track_idx = tracks_first[j].first
            particle_type = types[j]

            mass_calculated = INVALID_MASS  # Invalid marker

            # Handle cluster-based particles
            if cluster_idx < cluster_limit
                if particle_type == PDG_K_LONG  # Neutral hadron
                    nh_idx = cluster_idx + 1 - gamma_len
                    hit_idx = nhdata[nh_idx].hits.first + 1

                    # Get hit data
                    hit = calohits[hit_idx]
                    tof = hit.time

                    # Calculate distance
                    dx = hit.position.x - vx
                    dy = hit.position.y - vy
                    dz = hit.position.z - vz
                    l = sqrt(dx^2 + dy^2 + dz^2) * MM_TO_M

                    beta = l / (tof * C_LIGHT)

                    if 0.0f0 < beta < 1.0f0
                        energy = energies[j]
                        mass_calculated = energy * sqrt(1.0f0 - beta^2)
                    else
                        mass_calculated = INVALID_TOF_MASS  # Invalid measurement
                    end

                elseif particle_type == PDG_PHOTON  # Photon
                    mass_calculated = 0.0f0
                end
            end

            # Handle track-based particles (only if not already calculated)
            if mass_calculated < 0.0f0 && track_idx < tracks_len
                charge = charges[j]
                if abs(charge) > 0.0f0
                    mass = masses[j]

                    # Check for known particles
                    if abs(mass - ELECTRON_MASS) < ELECTRON_TOLERANCE
                        mass_calculated = ELECTRON_MASS
                    elseif abs(mass - MUON_MASS) < MUON_TOLERANCE
                        mass_calculated = MUON_MASS
                    else
                        # Calculate from time of flight
                        track = trackdata[track_idx + 1]
                        last_hit_idx = track.trackerHits.last
                        Tout = trackerhits[last_hit_idx].time
                        tof = Tout - v_t_scaled

                        l = track_L[track_idx + 1] * MM_TO_M
                        beta = l / (tof * C_LIGHT)

                        if 0.0f0 < beta < 1.0f0
                            # Calculate momentum magnitude
                            p = sqrt(mom_x[j]^2 + mom_y[j]^2 + mom_z[j]^2)
                            mass_calculated = p * sqrt(1.0f0 / (beta^2) - 1.0f0)
                        else
                            mass_calculated = PION_MASS  # Default
                        end
                    end
                end
            end

            # Store result if we calculated a mass
            if mass_calculated >= 0.0f0
                result_idx += 1
                tmp[result_idx] = mass_calculated
            end
        end

        # Resize to actual number of particles with calculated mass
        result[i] = resize!(tmp, result_idx)
    end

    return result
end

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

        @inbounds for j in 1:n_constituents
            d0_val = d0_vals[j]
            phi_val = phi_vals[j]

            sin_phi, cos_phi = sincos(phi_val)

            d0x = -d0_val * sin_phi
            d0y = d0_val * cos_phi

            dot_product = d0x * px + d0y * py

            abs_d0 = abs(d0_val)
            sign_dot = sign(dot_product)
            signed_ip = sign_dot * abs_d0

            is_valid = Float32(d0_val != UNDEF_VAL)
            sip2d_values[j] = is_valid * signed_ip + (1.0f0 - is_valid) * UNDEF_VAL
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

            sig_values[j] = valid ? sig : UNDEF_VAL
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

        @inbounds for j in 1:n_constituents
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

            is_valid = Float32(d0_val != UNDEF_VAL)
            cprojs[j] = is_valid * signed_ip + (1.0f0 - is_valid) * UNDEF_VAL
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

            sigs[j] = valid ? sig : UNDEF_VAL
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
                            jets_constituents::Vector{<:JetConstituents},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32) -> Vector{JetConstituentsData}

Calculate the jet distance value for each particle, measuring the distance between
the point of closest approach and the jet axis.
"""
function get_JetDistVal_clusterV(jets::Vector{JetReconstruction.EEJet},
                                 jets_constituents::Vector{<:JetConstituents},
                                 D0::Vector{JetConstituentsData},
                                 Z0::Vector{JetConstituentsData},
                                 phi0::Vector{JetConstituentsData},
                                 Bz::Float32)
    n_jets = length(jets)
    result = Vector{JetConstituentsData}(undef, n_jets)

    for i in 1:n_jets
        px_jet, py_jet, pz_jet = jets[i].px, jets[i].py, jets[i].pz
        single_jet_constituents = jets_constituents[i]
        n_constituents = length(single_jet_constituents)
        tmp = Vector{Float32}(undef, n_constituents)

        for j in 1:n_constituents
            if D0[i][j] != UNDEF_VAL
                d0_val = D0[i][j]
                z0_val = Z0[i][j]

                # Use sincos for efficiency
                sin_phi, cos_phi = sincos(phi0[i][j])

                # Impact parameter vector
                dx = -d0_val * sin_phi
                dy = d0_val * cos_phi
                dz = z0_val

                # Constituent momentum
                px_ct = single_jet_constituents[j].momentum.x
                py_ct = single_jet_constituents[j].momentum.y
                pz_ct = single_jet_constituents[j].momentum.z

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
                tmp[j] = UNDEF_VAL
            end
        end

        result[i] = tmp
    end

    return result
end

function get_btagJetDistVal(jets::Vector{JetReconstruction.EEJet},
                            jets_constituents::Vector{<:JetConstituents},
                            D0::Vector{JetConstituentsData},
                            Z0::Vector{JetConstituentsData},
                            phi0::Vector{JetConstituentsData},
                            Bz::Float32)
    # Simply call the implementation function
    return get_JetDistVal_clusterV(jets, jets_constituents, D0, Z0, phi0, Bz)
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

            tmp[j] = valid ? sig : UNDEF_VAL
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
