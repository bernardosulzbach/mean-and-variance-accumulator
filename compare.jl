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
mutable struct ReferenceAccumulator
    data::Vector{BigFloat}

    ReferenceAccumulator() = new(BigFloat[])
end

function add_sample!(acc::ReferenceAccumulator, x::Float32)
    push!(acc.data, BigFloat(x))
end

function get_mean(acc::ReferenceAccumulator)
    return Float64(sum(acc.data) / length(acc.data))
end

function get_variance(acc::ReferenceAccumulator)
    n = length(acc.data)
    mean_val = sum(acc.data) / n
    return Float64(sum((acc.data .- mean_val) .^ 2) / (n - 1))
end

function reset!(acc::ReferenceAccumulator)
    empty!(acc.data)
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
function display_results(name::String, mean_val, var_val, ref_mean, ref_var)
    println("$(name):")
    @printf("  Mean: %.10f, Variance: %.10f\n", mean_val, var_val)
    if ref_mean !== nothing && ref_var !== nothing
        @printf(
            "  Error - Mean: %+.2e, Variance: %+.2e\n",
            mean_val - ref_mean,
            var_val - ref_var
        )
    end
    println()
end

# Main comparison function
function compare_accumulators(data)
    # Create accumulator instances
    accumulators = [
        ("Naïve", NaiveAccumulator()),
        ("Welford (Welford, West, and Hanson)", WelfordAccumulatorWelfordWestHanson()),
        ("Welford (Youngs and Cramer)", WelfordAccumulatorYoungsCramer()),
        ("Welford (Knuth)", WelfordAccumulatorKnuth()),
    ]

    # High precision reference
    ref_acc = ReferenceAccumulator()
    ref_mean, ref_var = process_data!(ref_acc, data)

    println("Comparing variance accumulation methods")
    println("Data: $(length(data)) samples from Normal(1, 1)")
    println("Accumulator precision: Float32")
    println("=" ^ 80)

    display_results("Reference (high precision)", ref_mean, ref_var, nothing, nothing)

    # Process each accumulator
    results = []
    for (name, acc) in accumulators
        mean_val, var_val = process_data!(acc, data)
        display_results(name, mean_val, var_val, ref_mean, ref_var)
        push!(results, (name, mean_val, var_val))
    end

    # Summary of absolute errors
    println("Absolute errors:")
    for (name, mean_val, var_val) in results
        mean_err = abs(mean_val - ref_mean)
        var_err = abs(var_val - ref_var)
        @printf("%-40s - Mean: %.2e, Variance: %.2e\n", name, mean_err, var_err)
    end

    println("\nTheoretical values: Mean = 1.0, Variance = 1.0")
end

# Generate test data and run comparison
n_samples = 10 ^ 4
# Pick any specific generator for reproducibility
rng = Xoshiro(1)
# Normal(1, 1)
data = Float32.(randn(rng, n_samples) .+ 1.0)

compare_accumulators(data)
