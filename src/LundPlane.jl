"""
    decluster(jet::T, clusterseq::ClusterSequence{T}) where {T <: FourMomentum}

Given a jet and its clustering sequence, this function attempts to 
decluster it into its two parent subjets. If both parents exist, it 
returns them ordered by descending `pt²`.

Returns:
- `(j1, j2)` where `j1` is the harder subjet.
"""
function decluster(jet::T, clusterseq::ClusterSequence{T}) where {T <: FourMomentum}
    j1, j2 = parent_jets(jet, clusterseq)

    # Ensure both subjets exist
    if !(isnothing(j2) || isnothing(j1))
        # Order by descending pt²: j1 = harder, j2 = softer
        j1, j2 = pt2(j1) > pt2(j2) ? (j1, j2) : (j2, j1)
    end

    return (j1, j2)
end

"""
    generate_lund_emissions(jet::PseudoJet, cs::ClusterSequence{PseudoJet})

Generates the Lund plane emissions for a given jet. 
The jet is reclustered using the CA algorithm with a very
large R to fully capture the jet structure.

Returns:
- `lundPoints`: A vector of named tuples, each representing one step in the
  declustering with the following fields:
  - `h_pt`: harder branch pt
  - `s_pt`: softer branch pt
  - `z`: momentum fraction of the softer branch
  - `delta`: angular distance between branches
  - `kt`: transverse momentum of the softer branch relative to the harder
  - `psi`: azimuthal angle between branches
  - `kappa`: z * ΔR
"""
function generate_lund_emissions(jet::PseudoJet, cs::ClusterSequence{PseudoJet})

    # Recluster the input jet using Cambridge/Aachen with large R
    reconstructed_cluster_seq = recluster(jet, cs; algorithm = JetAlgorithm.CA,
                                          R = 1000.0)
    reconstructed_jet = inclusive_jets(reconstructed_cluster_seq; T = PseudoJet)[1]

    lundPoints = Vector{NamedTuple}()
    p1 = reconstructed_jet

    while true
        # Attempt to decluster p1 into p1 and p2
        p1, p2 = decluster(p1, reconstructed_cluster_seq)

        if isnothing(p1) || isnothing(p2)
            break  # No further parents
        end

        # Basic kinematic quantities
        harder_pt = pt(p1)
        softer_pt = pt(p2)

        Δ = deltaR(p1, p2)                            # Angular distance
        kt = pt(p2) * Δ                               # Relative transverse momentum
        z = pt(p2) / (pt(p1) + pt(p2))                # Momentum sharing
        psi = atan(rapidity(p2) - rapidity(p1), delta_phi(p2, p1))  # Azimuthal angle in (Δφ, Δy) plane
        kappa = z * Δ                                 # Log-polar projection

        # Store one step of the declustering
        temp = (h_pt = harder_pt,
                s_pt = softer_pt,
                z = z,
                delta = Δ,
                kt = kt,
                psi = psi,
                kappa = kappa)

        push!(lundPoints, temp)
    end

    return lundPoints
end
