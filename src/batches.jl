# The parallel modes split the paths into batches, each a serial run with its own seed and generator.
# Here is all they share: a mode only says how the batches run, with its own method of run_on_batches.

const BatchedConfig = MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractBatchedMode}

# A generator a batch can own. The extensions of generators that cannot be copied refuse here.
batch_rng(rng::Random.AbstractRNG) = copy(rng)

# Seeded here too: with init_rng = false the metrics do not, and the batches would share the same stream.
function batch_mode(rng::Random.AbstractRNG, seed::Int64)
    rng_ = batch_rng(rng)
    inner_seed!(rng_, seed)
    return SerialMode(seed, rng_)
end

function batch_modes(sub_mod::SerialMode, numberOfBatch::Int64, seeds::Array{Int64})
    ChainRulesCore.@ignore_derivatives @assert numberOfBatch > 0 "Number of batches must be positive"
    ChainRulesCore.@ignore_derivatives @assert length(seeds) == numberOfBatch "There must be one seed per batch"
    return [batch_mode(sub_mod.rng, seed) for seed in seeds]
end

paths_per_draw(::AbstractMonteCarloMethod) = 1
paths_per_draw(::AntitheticMC) = 2

# The batch sizes add up to Nsim: what does not divide evenly goes to the first batches
function batch_sizes(mcConfig::BatchedConfig)
    numberOfBatch = length(mcConfig.parallelMode.sub_modes)
    draws = paths_per_draw(mcConfig.monteCarloMethod)
    size_, remainder = divrem(div(mcConfig.Nsim, draws), numberOfBatch)
    ChainRulesCore.@ignore_derivatives size_ > 0 || throw(ArgumentError("$(mcConfig.Nsim) paths are not enough for $numberOfBatch batches"))
    return [draws * (size_ + (i <= remainder)) for i = 1:numberOfBatch]
end

# With init_rng = false the streams continue from a call to the next, in the generators of the batches
check_init_rng(::AbstractBatchedMode, ::Bool) = nothing
check_init_rng(::AbstractMultiProcess, init_rng::Bool) = ChainRulesCore.@ignore_derivatives init_rng || throw(ArgumentError("MultiProcess needs init_rng = true: the workers receive copies of the generators, so a stream cannot continue from a call to the next"))

function batch_configs(mcConfig::BatchedConfig)
    check_init_rng(mcConfig.parallelMode, mcConfig.init_rng)
    return [MonteCarloConfiguration(Nsim, mcConfig.Nstep, mcConfig.monteCarloMethod, mode, mcConfig.init_rng) for (Nsim, mode) in zip(batch_sizes(mcConfig), mcConfig.parallelMode.sub_modes)]
end

# The points of a Sobol sequence depend only on Nsim: the batches would all be the same
batch_configs(::MonteCarloConfiguration{<:Integer, <:Integer, <:SobolMode, <:AbstractBatchedMode}) = throw(ArgumentError("SobolMode does not support batches: every batch would simulate the same paths"))

function average_on_batches(f, mcConfig::BatchedConfig)
    configs = batch_configs(mcConfig)
    results = run_on_batches(f, configs, mcConfig.parallelMode)
    return sum(config.Nsim .* result for (config, result) in zip(configs, results)) ./ mcConfig.Nsim
end

set_seed!(mcConfig::BatchedConfig) = foreach(set_seed!, batch_configs(mcConfig))

# The paths of all the batches, in the order of the batches
join_paths(paths::Vector{<:AbstractMatrix}) = reduce(vcat, paths)
join_paths(paths::Vector{<:Vector{<:AbstractMatrix}}) = [reduce(vcat, [batch[i] for batch in paths]) for i in eachindex(paths[1])]

joined_on_batches(f, mcConfig::BatchedConfig) = join_paths(run_on_batches(f, batch_configs(mcConfig), mcConfig.parallelMode))

simulate(mcProcess::AbstractMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, T::Number) = joined_on_batches(config -> simulate(mcProcess, zeroCurve, config, T), mcConfig)
simulate(mcProcess::VectorialMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, T::Number) = joined_on_batches(config -> simulate(mcProcess, zeroCurve, config, T), mcConfig)
simulate(mcProcess::AbstractMonteCarloEngine, mcConfig::BatchedConfig, T::Number) = joined_on_batches(config -> simulate(mcProcess, config, T), mcConfig)

for (model_type, payoff_type) in ((BaseProcess, AbstractPayoff), (BaseProcess, Position), (BaseProcess, Array{<:AbstractPayoff}), (BaseProcess, AbstractSmile), (VectorialMonteCarloProcess, Array{Position}), (MarketDataSet, Portfolio))
    @eval function pricer(mcProcess::$model_type, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, abstractPayoff::$payoff_type)
        return average_on_batches(config -> pricer(mcProcess, rfCurve, config, abstractPayoff), mcConfig)
    end
end

# The main run of the control variates is what goes on the batches
function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:ControlVariates, <:AbstractBatchedMode}, abstractPayoff::AbstractPayoff)
    return pricer_control_variates(mcProcess, rfCurve, mcConfig, abstractPayoff)
end

# Mean and sum of the squared deviations of the discounted payoffs, one per payoff: what pools exactly over the batches
function path_moments(Payoff::AbstractVector)
    mean_ = mean(Payoff)
    return mean_, sum(abs2, Payoff .- mean_)
end

function path_moments(Payoffs::AbstractMatrix)
    means = mean(Payoffs, dims = 1)
    return vec(means), vec(sum(abs2, Payoffs .- means, dims = 1))
end

function batch_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::Union{AbstractPayoff, AbstractSmile})
    set_seed!(mcConfig)
    S = simulate(mcProcess, rfCurve, mcConfig, maturity(abstractPayoff))
    return path_moments(payoff(S, abstractPayoff, rfCurve, mcConfig))
end

function batch_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{<:AbstractPayoff})
    set_seed!(mcConfig)
    maxT = maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    moments = [path_moments(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]
    return first.(moments), last.(moments)
end

# The means and the variances of the paths of all the batches, without bringing the paths together
function pooled_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, abstractPayoffs)
    configs = batch_configs(mcConfig)
    moments = run_on_batches(config -> batch_moments(mcProcess, rfCurve, config, abstractPayoffs), configs, mcConfig.parallelMode)
    Means = sum(config.Nsim .* moment[1] for (config, moment) in zip(configs, moments)) ./ mcConfig.Nsim
    squares = sum(moment[2] .+ config.Nsim .* abs2.(moment[1] .- Means) for (config, moment) in zip(configs, moments))
    return Means, squares ./ (mcConfig.Nsim - 1)
end

for payoff_type in (AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile)
    @eval function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, abstractPayoffs::$payoff_type)
        return pooled_moments(mcProcess, rfCurve, mcConfig, abstractPayoffs)[2]
    end
    @eval function confinter(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, abstractPayoffs::$payoff_type, alpha::Real = 0.99)
        Means, Vars = pooled_moments(mcProcess, rfCurve, mcConfig, abstractPayoffs)
        IC = student_intervals(Means, Vars, mcConfig.Nsim, alpha)
        return Means isa Number ? only(IC) : IC
    end
end
