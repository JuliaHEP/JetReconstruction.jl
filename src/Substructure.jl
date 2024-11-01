# function to see if  a jet has parent jets
function has_parents(jet, history)
    
    N = jet._cluster_hist_index
    p1 = history[N].parent1
    p2 = history[N].parent2

    p1 == p2 == NonexistentParent ? result = false : result = true
    return result, p1, p2

end

# Function to calculate ΔR
function delta_R(jet1, jet2)
    
    eta1, phi1 = rapidity(jet1), phi(jet1)
    eta2, phi2 = rapidity(jet2), phi(jet2)

    d_eta = eta1 - eta2
    d_phi = phi1 - phi2
    d_phi = abs(d_phi) > π ? 2π - abs(d_phi) : d_phi
    
    return sqrt(d_eta^2 + d_phi^2)

end

# Function to calculate kt distance between two pseudojets
function kt_distance(jet1::PseudoJet, jet2::PseudoJet, R = 1)
    p1 = pt2(jet1)
    p2 = pt2(jet2)
    
    d_R = delta_R(jet1, jet2)
    return min(p1, p2) * (d_R^2 / R^2)
end;

# Function to recluster the constituents of a jet in a different method
function recluster(jet::PseudoJet, clusterseq::ClusterSequence, rad = 1.0, mtd = 0)
    cons = constituents(jet, clusterseq)
    new_clusterseq = jet_reconstruct(cons; p = mtd, R = rad, strategy = RecoStrategy.N2Plain)
   
    return new_clusterseq
end

function sort_jets!(event_jet_array)
    jet_pt(jet) = pt2(jet)
    sort!(event_jet_array, by = jet_pt, rev = true)
    # println(event_jet_array)
end

function join(jets::Vector{PseudoJet})
    final = PseudoJet(0.0, 0.0, 0.0, 0.0)
    
    for jet in jets
        final = final + jet
    end

    final
end

# Defining suitable structures to be used for jet grooming and tagging

struct MassDropTagger
    mu::Float64
    y::Float64
end

struct SoftDropTagger
    zcut::Float64
    b::Float64
end

struct Filter
    filterRadius::Float64 
    numHardestJets::Int
end

struct Trim
    trimRadius::Float64
    trimFraction::Float64
    reclusterMethod::Int
end

function massDrop(jet::PseudoJet, clusterseq::ClusterSequence, tag::MassDropTagger)
    allJets = clusterseq.jets
    hist = clusterseq.history

    while(true)
        had_parents, p1, p2 = has_parents(jet, hist)

        if (had_parents)
            parent1 = allJets[hist[p1].jetp_index]
            parent2 = allJets[hist[p2].jetp_index]

            if (m2(parent1) < m2(parent2))
                p1, p2 = p2, p1
                parent1, parent2 = parent2, parent1
            end

            if ((m2(parent1) < m2(jet)*tag.mu^2) && (kt_distance(parent1, parent2) > tag.y*m2(jet)))
                return jet
            else
                jet = parent1
            end
    
        else 
            return PseudoJet(0.0, 0.0, 0.0, 0.0)
        end

    end

end

function softDrop(jet::PseudoJet, clusterseq::ClusterSequence, rad::Float64, tag::SoftDropTagger)
    new_clusterseq = recluster(jet, clusterseq, rad, 0)
    new_jet = inclusive_jets(new_clusterseq; T = PseudoJet)[1]

    allJets = new_clusterseq.jets
    hist = new_clusterseq.history
    
    while(true)
        had_parents, p1, p2 = has_parents(new_jet, hist)

        if (had_parents)
            parent1 = allJets[hist[p1].jetp_index]
            parent2 = allJets[hist[p2].jetp_index]

            pti = pt2(parent1)^0.5
            ptj = pt2(parent2)^0.5

            if (m2(parent1) < m2(parent2))
                p1, p2 = p2, p1
                parent1, parent2 = parent2, parent1
            end

            if (min(pti, ptj)/(pti + ptj) > tag.zcut * (delta_R(parent1, parent2)/rad) ^ tag.b)
                return new_jet
            else
                new_jet = parent1
            end
    
        else 
            return PseudoJet(0.0, 0.0, 0.0, 0.0)
        end

    end

end

function jetFilter(jet::PseudoJet, clusterseq::ClusterSequence, filter::Filter)
    rad = filter.filterRadius;
    new_clusterseq = recluster(jet, clusterseq, rad, 0)
    reclustered = sort_jets!(inclusive_jets(new_clusterseq; T = PseudoJet))

    n = length(reclustered) <= filter.numHardestJets ? length(reclustered) : filter.numHardestJets 
    hard = reclustered[1:n]
    
    filtered = join(hard)
    
    filtered    
end

function jetTrim(jet::PseudoJet, clusterseq::ClusterSequence, trim::Trim)
    rad = trim.trimRadius;
    mtd = trim.reclusterMethod
    frac2 = trim.trimFraction ^ 2

    new_clusterseq = recluster(jet, clusterseq, rad, mtd)
    reclustered = sort_jets!(inclusive_jets(new_clusterseq; T = PseudoJet))
    
    hard = Vector{PseudoJet}(undef, 0)
    for item in reclustered
        if pt2(item) >= frac2 * pt2(jet)
            push!(hard, item)
        end
    end
    trimmed = join(hard)
    
    trimmed    
end;
