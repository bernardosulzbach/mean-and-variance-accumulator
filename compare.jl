using Random, Statistics, Printf

abstract type VarianceAccumulator end

mutable struct VarianceAccumulatorBase
    n::UInt32
    VarianceAccumulatorBase() = new(UInt32(0))
end

function add_sample!(acc::VarianceAccumulator, x::Float32) end
function get_mean(acc::VarianceAccumulator) end
function get_variance(acc::VarianceAccumulator) end
function reset!(acc::VarianceAccumulator) end

function increment_count!(acc::VarianceAccumulator)
    acc.base.n += UInt32(1)
end

function get_count(acc::VarianceAccumulator)
    return acc.base.n
end

function reset_count!(acc::VarianceAccumulator)
    acc.base.n = UInt32(0)
end

mutable struct NaiveAccumulator <: VarianceAccumulator
    base::VarianceAccumulatorBase
    sum_x::Float32
    sum_x2::Float32

    NaiveAccumulator() = new(VarianceAccumulatorBase(), Float32(0), Float32(0))
end

function add_sample!(acc::NaiveAccumulator, x::Float32)
    increment_count!(acc)

    acc.sum_x += x
    acc.sum_x2 += x * x
end

function get_mean(acc::NaiveAccumulator)
    return acc.sum_x / Float32(get_count(acc))
end

function get_variance(acc::NaiveAccumulator)
    n = Float32(get_count(acc))
    return (acc.sum_x2 - acc.sum_x * acc.sum_x / n) / (n - Float32(1))
end

function reset!(acc::NaiveAccumulator)
    reset_count!(acc)
    acc.sum_x = Float32(0)
    acc.sum_x2 = Float32(0)
end

mutable struct WelfordAccumulatorWelfordWestHanson <: VarianceAccumulator
    base::VarianceAccumulatorBase
    m::Float32
    s::Float32

    WelfordAccumulatorWelfordWestHanson() =
        new(VarianceAccumulatorBase(), Float32(0), Float32(0))
end

function add_sample!(acc::WelfordAccumulatorWelfordWestHanson, x::Float32)
    increment_count!(acc)

    j = Float32(get_count(acc))
    delta = x - acc.m
    acc.m += (1 / j) * delta
    acc.s += (j - 1) * delta * (delta / j)
end

function get_mean(acc::WelfordAccumulatorWelfordWestHanson)
    return acc.m
end

function get_variance(acc::WelfordAccumulatorWelfordWestHanson)
    n = Float32(get_count(acc))
    return acc.s / (n - Float32(1))
end

function reset!(acc::WelfordAccumulatorWelfordWestHanson)
    reset_count!(acc)
    acc.m = Float32(0)
    acc.s = Float32(0)
end

mutable struct WelfordAccumulatorWelfordWestHansonRearranged <: VarianceAccumulator
    base::VarianceAccumulatorBase
    m::Float32
    s::Float32

    WelfordAccumulatorWelfordWestHansonRearranged() =
        new(VarianceAccumulatorBase(), Float32(0), Float32(0))
end

function add_sample!(acc::WelfordAccumulatorWelfordWestHansonRearranged, x::Float32)
    increment_count!(acc)

    j = Float32(get_count(acc))
    delta = x - acc.m
    acc.m += delta / j
    acc.s += (j - 1) * delta * (delta / j)
end

function get_mean(acc::WelfordAccumulatorWelfordWestHansonRearranged)
    return acc.m
end

function get_variance(acc::WelfordAccumulatorWelfordWestHansonRearranged)
    n = Float32(get_count(acc))
    return acc.s / (n - Float32(1))
end

function reset!(acc::WelfordAccumulatorWelfordWestHansonRearranged)
    reset_count!(acc)
    acc.m = Float32(0)
    acc.s = Float32(0)
end

mutable struct WelfordAccumulatorYoungsCramer <: VarianceAccumulator
    base::VarianceAccumulatorBase
    t::Float32
    s::Float32

    WelfordAccumulatorYoungsCramer() =
        new(VarianceAccumulatorBase(), Float32(0), Float32(0))
end

function add_sample!(acc::WelfordAccumulatorYoungsCramer, x::Float32)
    increment_count!(acc)

    j = Float32(get_count(acc))
    acc.t += x
    if j > 1
        acc.s += (1 / (j * (j - 1))) * (j * x - acc.t) ^ 2
    end
end

function get_mean(acc::WelfordAccumulatorYoungsCramer)
    return acc.t / Float32(get_count(acc))
end

function get_variance(acc::WelfordAccumulatorYoungsCramer)
    n = Float32(get_count(acc))
    return acc.s / (n - Float32(1))
end

function reset!(acc::WelfordAccumulatorYoungsCramer)
    reset_count!(acc)
    acc.t = Float32(0)
    acc.s = Float32(0)
end

mutable struct WelfordAccumulatorYoungsCramerRearranged <: VarianceAccumulator
    base::VarianceAccumulatorBase
    t::Float32
    s::Float32

    WelfordAccumulatorYoungsCramerRearranged() =
        new(VarianceAccumulatorBase(), Float32(0), Float32(0))
end

function add_sample!(acc::WelfordAccumulatorYoungsCramerRearranged, x::Float32)
    increment_count!(acc)

    j = Float32(get_count(acc))
    acc.t += x
    if j > 1
        acc.s += (j * x - acc.t) ^ 2 / (j * (j - 1))
    end
end

function get_mean(acc::WelfordAccumulatorYoungsCramerRearranged)
    return acc.t / Float32(get_count(acc))
end

function get_variance(acc::WelfordAccumulatorYoungsCramerRearranged)
    n = Float32(get_count(acc))
    return acc.s / (n - Float32(1))
end

function reset!(acc::WelfordAccumulatorYoungsCramerRearranged)
    reset_count!(acc)
    acc.t = Float32(0)
    acc.s = Float32(0)
end

mutable struct WelfordAccumulatorKnuth <: VarianceAccumulator
    base::VarianceAccumulatorBase
    mean::Float32
    m2::Float32

    WelfordAccumulatorKnuth() = new(VarianceAccumulatorBase(), Float32(0), Float32(0))
end

function add_sample!(acc::WelfordAccumulatorKnuth, x::Float32)
    increment_count!(acc)

    n = Float32(get_count(acc))
    delta = x - acc.mean
    acc.mean += delta / n
    delta2 = x - acc.mean
    acc.m2 += delta * delta2
end

function get_mean(acc::WelfordAccumulatorKnuth)
    return acc.mean
end

function get_variance(acc::WelfordAccumulatorKnuth)
    n = Float32(get_count(acc))
    return acc.m2 / (n - Float32(1))
end

function reset!(acc::WelfordAccumulatorKnuth)
    reset_count!(acc)
    acc.mean = Float32(0)
    acc.m2 = Float32(0)
end

# High precision reference accumulator
mutable struct ReferenceAccumulator <: VarianceAccumulator
    base::VarianceAccumulatorBase
    sum_x::BigFloat
    sum_x2::BigFloat

    ReferenceAccumulator() = new(VarianceAccumulatorBase(), BigFloat(0), BigFloat(0))
end

function add_sample!(acc::ReferenceAccumulator, x::Float32)
    increment_count!(acc)

    acc.sum_x += BigFloat(x)
    acc.sum_x2 += BigFloat(x) * BigFloat(x)
end

function get_mean(acc::ReferenceAccumulator)
    return acc.sum_x / BigFloat(get_count(acc))
end

function get_variance(acc::ReferenceAccumulator)
    n = BigFloat(get_count(acc))
    return (acc.sum_x2 - acc.sum_x * acc.sum_x / n) / (n - BigFloat(1))
end

function reset!(acc::ReferenceAccumulator)
    reset_count!(acc)
    acc.sum_x = BigFloat(0)
    acc.sum_x2 = BigFloat(0)
end

# Generic function to process data with any accumulator
function process_data!(acc, data)
    reset!(acc)
    for x in data
        add_sample!(acc, x)
    end
    return get_mean(acc), get_variance(acc)
end

# Results display helper
function display_results(i, name::String, mean_val, var_val, ref_mean, ref_var)
    println("$(i). $(name):")
    @printf("   Mean: %.10f, Variance: %.10f\n", mean_val, var_val)
    if ref_mean !== nothing && ref_var !== nothing
        @printf(
            "   Error - Mean: %+.2e, Variance: %+.2e\n",
            mean_val - ref_mean,
            var_val - ref_var
        )
    end
    println()
end

# Main comparison function
function compare_accumulators()
    # Create accumulator instances
    accumulators = [
        ("Naïve", NaiveAccumulator()),
        ("Welford (Welford, West, and Hanson)", WelfordAccumulatorWelfordWestHanson()),
        (
            "Welford (Welford, West, and Hanson) rearranged",
            WelfordAccumulatorWelfordWestHansonRearranged(),
        ),
        ("Welford (Youngs and Cramer)", WelfordAccumulatorYoungsCramer()),
        (
            "Welford (Youngs and Cramer) rearranged",
            WelfordAccumulatorYoungsCramerRearranged(),
        ),
        ("Welford (Knuth)", WelfordAccumulatorKnuth()),
    ]

    ref_acc = ReferenceAccumulator()

    sample_counts = [10 ^ e for e = 1:9]
    function iterations_for(sample_count)
        return max(10, maximum(sample_counts) ÷ sample_count)
    end
    means = [
        zeros(Float64, iterations_for(sample_count)) for
        _ = 1:length(accumulators), sample_count in sample_counts
    ]
    vars = [
        zeros(Float64, iterations_for(sample_count)) for
        _ = 1:length(accumulators), sample_count in sample_counts
    ]
    for (i_sample_counts, sample_count) in enumerate(sample_counts)
        for i_iteration = 1:iterations_for(sample_count)
            # Use a specific PRNG to ensure reproducibility
            rng = Xoshiro(i_iteration * sample_count)
            # Normal(0.1, 0.01)
            data = randn(rng, Float32, sample_count) .* Float32(0.01) .+ Float32(0.1)
            ref_mean, ref_var = process_data!(ref_acc, data)
            for i_acc = 1:length(accumulators)
                _, acc = accumulators[i_acc]
                mean_val, var_val = process_data!(acc, data)
                means[i_acc, i_sample_counts][i_iteration] = abs(mean_val - ref_mean)
                vars[i_acc, i_sample_counts][i_iteration] = abs(var_val - ref_var)
            end
        end
    end

    for (i_sample_counts, sample_count) in enumerate(sample_counts)
        println("=" ^ 164)
        @printf(
            "%-48s | %-55s | %-55s\n",
            @sprintf(
                "%d samples (%d iterations)",
                sample_count,
                iterations_for(sample_count)
            ),
            "Absolute error of the mean",
            "Absolute error of the variance"
        )
        @printf(
            "%-48s | %-26s | %-26s | %-26s | %-26s\n",
            "Algorithm",
            "Mean (/ max)",
            "s (/ max)",
            "Mean (/ max)",
            "s (/ max)"
        )
        println("-" ^ 164)
        means_m = [mean(means[i, i_sample_counts]) for i = 1:length(accumulators)]
        means_s = [std(means[i, i_sample_counts]) for i = 1:length(accumulators)]
        vars_m = [mean(vars[i, i_sample_counts]) for i = 1:length(accumulators)]
        vars_s = [std(vars[i, i_sample_counts]) for i = 1:length(accumulators)]
        for (i, (name, _)) in enumerate(accumulators)
            @printf(
                "%-48s | %.9e (%.6f) | %.9e (%.6f) | %.9e (%.6f) | %.9e (%.6f)\n",
                name,
                means_m[i],
                means_m[i] / maximum(means_m),
                means_s[i],
                means_s[i] / maximum(means_s),
                vars_m[i],
                vars_m[i] / maximum(vars_m),
                vars_s[i],
                vars_s[i] / maximum(vars_s)
            )
        end
    end
end

compare_accumulators()
