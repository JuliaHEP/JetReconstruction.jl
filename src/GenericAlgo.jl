# This is the generic reconstruction algorithm that will
# switch based on the strategy, or based on the event density
# if the "Best" strategy is to be employed

function jet_reconstruct(particles; p = -1, R = 1.0, recombine = +, strategy = RecoStrategy.Best)
    # Either map to the fixed algorithm corresponding to the strategy
    # or to an optimal choice based on the density of initial particles

    if strategy == RecoStrategy.Best
        # The breakpoint of ~90 is determined empirically on e+e- -> H and 0.5TeV pp -> 5GeV jets
        algorithm = length(particles) > 80 ? tiled_jet_reconstruct : plain_jet_reconstruct
    elseif strategy == RecoStrategy.N2Plain
        algorithm = plain_jet_reconstruct
    elseif strategy == RecoStrategy.N2Tiled
        algorithm = tiled_jet_reconstruct
    else
        throw(ErrorException("Invalid strategy: $(strategy)"))
    end

    # Now call the chosen algorithm, passing through the other parameters
    algorithm(particles; p = p, R = R, recombine = recombine)
end
