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
    generate_lund_projection(jet::PseudoJet, cs::ClusterSequence{PseudoJet})

Generates the Lund plane projection for a given jet. 
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
function generate_lund_projection(jet::PseudoJet, cs::ClusterSequence{PseudoJet})

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

"""
generate_average_lund_image(njets::Int, delta_array::Vector{Vector{Real}}, 
                            kt_array::Vector{Vector{Real}}; xrange::Tuple{Real}, 
                            yrange::Tuple{Real}, bins::Int) -> (xgrid, ygrid, avg_image)

Computes an average Lund image from a set of jets by binning the (log(1/ΔR), log(kt))
coordinates into a fixed-size 2D histogram. Each jet's histogram is normalized and
then averaged across all jets.

Arguments:
- `njets`: Number of jets
- `delta_array`: Vector of vectors, where each inner vector contains ΔR values for a jet
- `kt_array`: Vector of vectors, where each inner vector contains kt values for a jet
- `xrange`: Tuple defining the x-axis (log(1/ΔR)) range
- `yrange`: Tuple defining the y-axis (log(kt)) range
- `bins`: Number of bins along each axis

Returns:
- A tuple `(xgrid, ygrid, avg_image)` where `xgrid` and `ygrid` are coordinate axes labels and
  `avg_image` is the averaged 2D histogram as a matrix
"""
function generate_average_lund_image(njets::Int, delta_array::Vector{Vector{Real}},
                                     kt_array::Vector{Vector{Real}}; xrange = (0.0, 9.0),
                                     yrange = (-5.0, 7.0), bins = 25)
    xmin, xmax = xrange
    ymin, ymax = yrange

    x_width = (xmax - xmin)/bins
    y_width = (ymax - ymin)/bins

    total_res = []

    for i in 1:njets
        Xind = ceil.((delta_array[i] .- xmin) ./ x_width)
        Yind = ceil.((kt_array[i] .- ymin) ./ y_width)

        res = zeros(Float32, bins, bins)
        L1norm = 0.0

        for i in 1:length(Xind)
            x = Int(Xind[i])
            y = Int(Yind[i])
            if (maximum([x, y]) < bins && minimum([x, y]) >= 1)
                res[x, y] += 1
                L1norm += 1
            end
        end

        if L1norm > 0
            res[1, :] .= res[1, :] ./ L1norm
            push!(total_res, res)
        end
    end

    avg_res = mean(total_res, dims = 1)[1]
    x = range(xmin, xmax; length = bins)
    y = range(ymin, ymax; length = bins)

    return (x, y, avg_res)
end
