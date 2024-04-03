# This is the generic reconstruction algorithm that will
# switch based on the strategy, or based on the event density
# if the "Best" strategy is to be employed

function generic_jet_reconstruct(particles; p = -1, R = 1.0, recombine = +, ptmin = 0.0, strategy = JetRecoStrategy.Best)
    # Either map to the fixed algorithm corresponding to the strategy
    # or to an optimal choice based on the density of initial particles

    if strategy == JetRecoStrategy.Best
        # The breakpoint of ~90 is determined empirically on e+e- -> H and 0.5TeV pp -> 5GeV jets
        algorithm = length(particles) > 80 ? tiled_jet_reconstruct : plain_jet_reconstruct
    elseif strategy == JetRecoStrategy.N2Plain
        algorithm = plain_jet_reconstruct
    elseif strategy == JetRecoStrategy.N2Tiled
        algorithm = tiled_jet_reconstruct
    else
        throw(ErrorException("Invalid strategy: $(strategy)"))
    end

    # Now call the chosen algorithm, passing through the other parameters
    algorithm(particles; p = p, R = R, recombine = recombine, ptmin = ptmin)
end
