using EDM4hep
using JetReconstruction
using StructArrays: StructVector

const C_LIGHT = 2.99792458e8
const C_LIGHT_INV = 1.0f0 / C_LIGHT
const ELECTRON_MASS = 0.000510999f0
const MUON_MASS = 0.105658f0
const PION_MASS = 0.13957039f0
const ELECTRON_TOLERANCE = 1.0f-5
const MUON_TOLERANCE = 1.0f-3

### Special Measurements (2)

# get_dndx - dE/dx measurement (energy loss)
# get_mtof - Mass from time-of-flight measurement

"""
    get_dndx(jcs::Vector{JetConstituents}, 
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
function get_dndx(jcs::Vector{<:JetConstituents}, 
                  dNdx::StructVector{EDM4hep.Quantity},
                  trackdata::StructVector{EDM4hep.Track}, 
                  JetsConstituents_isChargedHad::Vector{Vector{Float32}})
    
    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    tracks_len = length(trackdata)

    @inbounds for i in 1:n_jets
        jet = jcs[i]
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

function get_mtof(jcs::Vector{<:JetConstituents}, 
                  track_L::AbstractArray{T} where T <: Float32,
                  trackdata::StructVector{EDM4hep.Track},
                  trackerhits::StructVector{EDM4hep.TrackerHit},
                  gammadata::StructVector{EDM4hep.Cluster},
                  nhdata::StructVector{EDM4hep.Cluster},
                  calohits::StructVector{EDM4hep.CalorimeterHit},
                  V::LorentzVector)
    
    n_jets = length(jcs)
    result = Vector{Vector{Float32}}(undef, n_jets)
    
    # Pre-compute limits
    tracks_len = length(trackdata)
    gamma_len = length(gammadata)
    nh_len = length(nhdata)
    cluster_limit = nh_len + gamma_len
    
    # Pre-compute vertex values
    vx, vy, vz = V.x, V.y, V.z
    v_t_scaled = V.t * 1.0f-3 * C_LIGHT_INV  # Tin calculation
    
    @inbounds for i in 1:n_jets
        jc = jcs[i]
        n_constituents = length(jc)
        
        # Pre-allocate for this jet
        tmp = Vector{Float32}(undef, n_constituents)
        result_idx = 0
        
        # Access fields once
        clusters_first = jc.clusters
        tracks_first = jc.tracks
        types = jc.type
        charges = jc.charge
        masses = jc.mass
        energies = jc.energy
        mom_x = jc.momentum.x
        mom_y = jc.momentum.y
        mom_z = jc.momentum.z
        
        for j in 1:n_constituents
            cluster_idx = clusters_first[j].first
            track_idx = tracks_first[j].first
            particle_type = types[j]
            
            mass_calculated = -1.0f0  # Invalid marker
            
            # Handle cluster-based particles
            if cluster_idx < cluster_limit
                if particle_type == 130  # Neutral hadron
                    nh_idx = cluster_idx + 1 - gamma_len
                    hit_idx = nhdata[nh_idx].hits.first + 1
                    
                    # Get hit data
                    hit = calohits[hit_idx]
                    tof = hit.time
                    
                    # Calculate distance
                    dx = hit.position.x - vx
                    dy = hit.position.y - vy
                    dz = hit.position.z - vz
                    L = sqrt(dx^2 + dy^2 + dz^2) * 0.001f0
                    
                    beta = L / (tof * C_LIGHT)
                    
                    if 0.0f0 < beta < 1.0f0
                        E = energies[j]
                        mass_calculated = E * sqrt(1.0f0 - beta^2)
                    else
                        mass_calculated = 9.0f0  # Invalid measurement
                    end
                    
                elseif particle_type == 22  # Photon
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
                        
                        L = track_L[track_idx + 1] * 0.001f0
                        beta = L / (tof * C_LIGHT)
                        
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