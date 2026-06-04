using Statistics: mean, std, median, quantile
using Printf: @sprintf


function summarise_trial_times(trial_timing::AbstractVector{<:Real})
    (
        n_samples = length(trial_timing),
        mean = mean(trial_timing),
        std = length(trial_timing) > 1 ? std(trial_timing) : 0.0,
        median = median(trial_timing),
        minimum = minimum(trial_timing),
        maximum = maximum(trial_timing),
        q25 = quantile(trial_timing, 0.25),
        q75 = quantile(trial_timing, 0.75),
        iqr = quantile(trial_timing, 0.75) - quantile(trial_timing, 0.25),
    )
end

function filter_outliers_iqr(trial_timing::AbstractVector{<:Real}, stats; outlier_band::Real = 2.0)
    min_val = stats.q25 - outlier_band * stats.iqr
    max_val = stats.q75 + outlier_band * stats.iqr
    trial_timing[(time -> min_val <= time <= max_val).(trial_timing)]
end

function pprint_trial_stats(stats)
    " - average time per event " * @sprintf("%.2f", stats.mean) * " ± " *
    @sprintf("%.2f", stats.std) * " μs\n" *
    " - median time per event " * @sprintf("%.2f", stats.median) * " μs\n" *
    " - lowest time per event " * @sprintf("%.2f", stats.minimum) * " μs"
end

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
