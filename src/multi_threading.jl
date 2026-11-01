abstract type AbstractMultiThreading <: ParallelMode end

# A generator a batch can own. The extensions of generators that cannot be copied refuse here.
batch_rng(rng::Random.AbstractRNG) = copy(rng)

# Seeded here too: with init_rng = false the metrics do not, and the batches would share the same stream.
function batch_mode(rng::Random.AbstractRNG, seed::Int64)
    rng_ = batch_rng(rng)
    inner_seed!(rng_, seed)
    return SerialMode(seed, rng_)
end

"""
    MultiThreading(sub_mod = SerialMode(), numberOfBatch = Threads.nthreads(), seeds = sub_mod.seed .+ (1:numberOfBatch))

Parallel mode splitting the paths into `numberOfBatch` batches priced on different threads, each
with its own seed from `seeds` and its own copy of the generator of `sub_mod`, a [`FinancialMonteCarlo.SerialMode`](@ref).
The result does not depend on the number of threads. The price is the average of the batch prices,
weighted by their number of paths, and the other metrics work on the paths of all the batches.

Generators that cannot be copied (VectorizedRNG) and [`FinancialMonteCarlo.SobolMode`](@ref) are not
supported. Of the automatic differentiation backends, only the forward ones are supported: the
batches would write to the same tape of ReverseDiff, and Mooncake cannot differentiate through tasks.
"""
struct MultiThreading{modeType <: SerialMode} <: AbstractMultiThreading
    sub_modes::Vector{modeType}
    function MultiThreading(sub_mod::SerialMode = SerialMode(), numberOfBatch::Int64 = Int64(Threads.nthreads()), seeds::Array{Int64} = sub_mod.seed .+ collect(1:numberOfBatch))
        ChainRulesCore.@ignore_derivatives @assert numberOfBatch > 0 "Number of batches must be positive"
        ChainRulesCore.@ignore_derivatives @assert length(seeds) == numberOfBatch "There must be one seed per batch"
        sub_modes = [batch_mode(sub_mod.rng, seed) for seed in seeds]
        return new{eltype(sub_modes)}(sub_modes)
    end
end

const MultiThreadingConfig = MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractMultiThreading}

paths_per_draw(::AbstractMonteCarloMethod) = 1
paths_per_draw(::AntitheticMC) = 2

# The batch sizes add up to Nsim: what does not divide evenly goes to the first batches
function batch_sizes(mcConfig::MultiThreadingConfig)
    numberOfBatch = length(mcConfig.parallelMode.sub_modes)
    draws = paths_per_draw(mcConfig.monteCarloMethod)
    size_, remainder = divrem(div(mcConfig.Nsim, draws), numberOfBatch)
    ChainRulesCore.@ignore_derivatives size_ > 0 || throw(ArgumentError("$(mcConfig.Nsim) paths are not enough for $numberOfBatch batches"))
    return [draws * (size_ + (i <= remainder)) for i = 1:numberOfBatch]
end

function batch_configs(mcConfig::MultiThreadingConfig)
    return [MonteCarloConfiguration(Nsim, mcConfig.Nstep, mcConfig.monteCarloMethod, mode, mcConfig.init_rng) for (Nsim, mode) in zip(batch_sizes(mcConfig), mcConfig.parallelMode.sub_modes)]
end

# The points of a Sobol sequence depend only on Nsim: the batches would all be the same
batch_configs(::MonteCarloConfiguration{<:Integer, <:Integer, <:SobolMode, <:AbstractMultiThreading}) = throw(ArgumentError("SobolMode does not support MultiThreading: every batch would simulate the same paths"))

run_on_batches(f, configs) = fetch.([Threads.@spawn f(config) for config in configs])

function average_on_batches(f, mcConfig::MultiThreadingConfig)
    configs = batch_configs(mcConfig)
    results = run_on_batches(f, configs)
    return sum(config.Nsim .* result for (config, result) in zip(configs, results)) ./ mcConfig.Nsim
end

set_seed!(mcConfig::MultiThreadingConfig) = foreach(set_seed!, batch_configs(mcConfig))

# The paths of all the batches, in the order of the batches: what the metrics not defined here work on
join_paths(paths::Vector{<:AbstractMatrix}) = reduce(vcat, paths)
join_paths(paths::Vector{<:Vector{<:AbstractMatrix}}) = [reduce(vcat, [batch[i] for batch in paths]) for i in eachindex(paths[1])]

simulate(mcProcess::AbstractMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcConfig::MultiThreadingConfig, T::Number) = join_paths(run_on_batches(config -> simulate(mcProcess, zeroCurve, config, T), batch_configs(mcConfig)))
simulate(mcProcess::VectorialMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcConfig::MultiThreadingConfig, T::Number) = join_paths(run_on_batches(config -> simulate(mcProcess, zeroCurve, config, T), batch_configs(mcConfig)))
simulate(mcProcess::AbstractMonteCarloEngine, mcConfig::MultiThreadingConfig, T::Number) = join_paths(run_on_batches(config -> simulate(mcProcess, config, T), batch_configs(mcConfig)))

for (model_type, payoff_type) in ((BaseProcess, AbstractPayoff), (BaseProcess, Position), (BaseProcess, Array{<:AbstractPayoff}), (BaseProcess, AbstractSmile), (VectorialMonteCarloProcess, Array{Position}), (MarketDataSet, Portfolio))
    @eval function pricer(mcProcess::$model_type, rfCurve::AbstractZeroRateCurve, mcConfig::MultiThreadingConfig, abstractPayoff::$payoff_type)
        return average_on_batches(config -> pricer(mcProcess, rfCurve, config, abstractPayoff), mcConfig)
    end
end

# The main run of the control variates is what goes on the batches
function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:ControlVariates, <:AbstractMultiThreading}, abstractPayoff::AbstractPayoff)
    return pricer_control_variates(mcProcess, rfCurve, mcConfig, abstractPayoff)
end
