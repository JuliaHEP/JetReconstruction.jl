module JetBenchmarkPlots

using JetReconstruction
using UnicodePlots

"""
    plot_trial_times(trial_timing)

Plot terminal diagnostics for per-event benchmark timing measurements.

This method is provided by the `UnicodePlots` extension. It prints a histogram
of trial timings and a line plot showing how the timing varies across trials.
"""
function JetReconstruction.plot_trial_times(trial_timing::AbstractVector{<:Real})
    stats = JetReconstruction.summarise_trial_times(trial_timing)

    # Freedman-Diaconis bin width, using IQR to remain robust in the presence of
    # outlier timings from system noise.
    bin_width = ceil(2 * stats.iqr / stats.n_samples^(1 / 3))

    if bin_width <= 0 || !isfinite(bin_width)
        bin_width = 1
    end

    nbins = min(max(ceil(Int, (stats.maximum - stats.minimum) / bin_width), 1), 80)

    println(
        UnicodePlots.histogram(
            trial_timing;
            nbins = nbins,
            vertical = true,
            title = "Histogram of event time per trial",
        ),
    )

    println(
        UnicodePlots.lineplot(
            collect(1:length(trial_timing)),
            trial_timing;
            title = "Runtime per event across trials",
        ),
    )

    return nothing
end

end
