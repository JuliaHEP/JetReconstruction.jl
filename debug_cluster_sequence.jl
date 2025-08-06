# Print cluster sequence and jet info for comparison with FastJet
using Printf
import JetReconstruction: pt, rapidity, px, py, pz, energy

function debug_cluster_sequence(clusterseq, particles, label="ClusterSeq")
    println("=== CLUSTER SEQUENCE DEBUG ===")
    println("Initial particles: ", length(particles))
    for (i, p) in enumerate(particles)
        pj = typeof(p) <: NamedTuple ? p : p
        println(@sprintf("Initial[%d]: px=%.6f, py=%.6f, pz=%.6f, E=%.6f, pt=%.6f, rap=%.7f", i-1, px(pj), py(pj), pz(pj), energy(pj), pt(pj), rapidity(pj)))
    end
    println()
    println("Cluster sequence history:")
    println("Total jets in sequence: ", length(clusterseq.jets))
    for (i, jet) in enumerate(clusterseq.jets)
        pj = typeof(jet) <: NamedTuple ? jet : jet
        origin = i <= length(particles) ? "(original particle)" : "(merged jet)"
        println(@sprintf("Jet[%d]: px=%.6f, py=%.6f, pz=%.6f, E=%.6f, pt=%.6f, rap=%.7f %s", i-1, px(pj), py(pj), pz(pj), energy(pj), pt(pj), rapidity(pj), origin))
    end
    println()
    if hasproperty(clusterseq, :history)
        println("Cluster sequence history steps:")
        println("History size: ", length(clusterseq.history))
        for (i, h) in enumerate(clusterseq.history)
            if hasproperty(h, :type) && h.type == :initial
                j = h.jet
                pj = clusterseq.jets[j+1]
                println(@sprintf("Step[%d]: INITIAL - particle %d, dij=0 (result: pt=%.6f, rap=%.7f)", i-1, j, pt(pj), rapidity(pj)))
            elseif hasproperty(h, :type) && h.type == :merge
                j1, j2, jnew = h.j1, h.j2, h.jnew
                dij = h.dij
                pj = clusterseq.jets[jnew+1]
                println(@sprintf("Step[%d]: JET merge - jets %d + %d -> %d, dij=%.6f (result: pt=%.6f, rap=%.7f)", i-1, j1, j2, jnew, dij, pt(pj), rapidity(pj)))
            elseif hasproperty(h, :type) && h.type == :beam
                j = h.jet
                dij = h.dij
                println(@sprintf("Step[%d]: BEAM merge - jet %d -> final, dij=%.6f", i-1, j, dij))
            end
        end
    end
    println()
    println("=== END CLUSTER SEQUENCE DEBUG ===")
end
