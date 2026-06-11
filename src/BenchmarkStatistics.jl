using Statistics: mean, std, median, quantile
using Printf: @sprintf

"""
    summarise_trial_times(trial_timing)

Summarise per-event timing measurements from repeated benchmark trials.

Returns a named tuple containing the number of samples, mean, standard
deviation, median, minimum, maximum, quartiles, and interquartile range.
"""
function summarise_trial_times(trial_timing::AbstractVector{<:Real})
    (n_samples = length(trial_timing),
     mean = mean(trial_timing),
     std = length(trial_timing) > 1 ? std(trial_timing) : 0.0,
     median = median(trial_timing),
     minimum = minimum(trial_timing),
     maximum = maximum(trial_timing),
     q25 = quantile(trial_timing, 0.25),
     q75 = quantile(trial_timing, 0.75),
     iqr = quantile(trial_timing, 0.75) - quantile(trial_timing, 0.25))
end

"""
    filter_outliers_iqr(trial_timing, stats; outlier_band=2.0)

Return the entries in `trial_timing` that lie within `outlier_band` times the
interquartile range around the first and third quartiles in `stats`.
"""
function filter_outliers_iqr(trial_timing::AbstractVector{<:Real}, stats;
                             outlier_band::Real = 2.0)
    min_val = stats.q25 - outlier_band * stats.iqr
    max_val = stats.q75 + outlier_band * stats.iqr
    trial_timing[(time -> min_val <= time <= max_val).(trial_timing)]
end

"""
    pprint_trial_stats(stats)

Format a timing summary, as returned by [`summarise_trial_times`](@ref), for
human-readable terminal output.
"""
function pprint_trial_stats(stats)
    " - average time per event " * @sprintf("%.2f", stats.mean) * " ± " *
    @sprintf("%.2f", stats.std) * " μs\n" *
    " - median time per event " * @sprintf("%.2f", stats.median) * " μs\n" *
    " - lowest time per event " * @sprintf("%.2f", stats.minimum) * " μs"
end

"""
    print_statistics(trial_timing; outlier_exclusion=true, outlier_band=2.0)

Print summary statistics for per-event benchmark timing measurements.

When `outlier_exclusion` is `true`, also print a second summary after excluding
measurements outside `outlier_band` times the interquartile range. Returns a
named tuple with the full summary and, when requested, the outlier-excluded
summary.
"""
function print_statistics(trial_timing::AbstractVector{<:Real};
                          outlier_exclusion::Bool = true,
                          outlier_band::Real = 2.0)
    stats = summarise_trial_times(trial_timing)
    println("Full statistics ($(stats.n_samples) samples)")
    println(pprint_trial_stats(stats))

    filtered_stats = nothing
    if outlier_exclusion
        no_outliers = filter_outliers_iqr(trial_timing, stats; outlier_band = outlier_band)
        filtered_stats = summarise_trial_times(no_outliers)
        println("Excluding outliers at $(outlier_band)xIQR (leaving $(filtered_stats.n_samples) of $(stats.n_samples) samples)")
        println(pprint_trial_stats(filtered_stats))
    end

    (full = stats, outlier_excluded = filtered_stats)
end
