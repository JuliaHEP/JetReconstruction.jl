# This is the generic reconstruction algorithm that will
# switch based on the strategy, or based on the event density
# if the "Best" strategy is to be employed

function generic_jet_reconstruct(particles; p = -1, R = 1.0, recombine = +, ptmin = 0.0, strategy = Best::JetRecoStrategy)
    # Either map to the fixed algorithm corresponding to the strategy
    # or to an optimal choice based on the density of initial particles
    algorithm = @match strategy begin
        N2Plain => plain_jet_reconstruct
        N2Tiled => tiled_jet_reconstruct
        Best => length(particles) > 50 ? tiled_jet_reconstruct : plain_jet_reconstruct
    end

    # Now call the chosen algorithm, passing through the other parameters
    algorithm(particles; p = p, R = R, recombine = recombine, ptmin = ptmin)
end
