# The parallel modes split the paths into batches, each a serial run with its own seed and generator.
# Here is all they share: a mode only says how the batches run, with its own method of run_on_batches.

paths_per_draw(::AbstractSampling) = 1
paths_per_draw(::AntitheticMC) = 2

# the message is built only when the check fails
check_batch_size(size_, Nsim, nb) = size_ > 0 || throw_argument_error("$Nsim paths are not enough for $nb batches")
ChainRulesCore.@non_differentiable check_batch_size(::Any, ::Any, ::Any)

# The batch sizes add up to Nsim: what does not divide evenly goes to the first batches
function batch_sizes(mcConfig::BatchedConfig)
    numberOfBatch = nbatches(mcConfig.mode)
    draws = paths_per_draw(mcConfig.sampling)
    size_, remainder = divrem(div(mcConfig.Nsim, draws), numberOfBatch)
    check_batch_size(size_, mcConfig.Nsim, numberOfBatch)
    return [draws * (size_ + (i <= remainder)) for i = 1:numberOfBatch]
end

# Each batch is a serial configuration owning the generator of its batch, so that with reseed = false the streams
# continue from a call to the next.
batch_configs(mcConfig::BatchedConfig) = [MonteCarloConfiguration(Nsim, mcConfig.Nstep, mcConfig.sampling, SerialMode(), state, mcConfig.reseed) for (Nsim, state) in zip(batch_sizes(mcConfig), mcConfig.rng)]

# The results of run_on_batches, with the type inference gives to f: fetch and pmap return Any, which would make every
# metric on batches type unstable
function batch_results(f, configs, mode::AbstractBatchedMode)
    results = run_on_batches(f, configs, mode)
    result_type = Base.promote_op(f, eltype(configs))
    return result_type === Union{} ? results : convert(Vector{result_type}, results)
end

weighted_average(configs, results, mcConfig::BatchedConfig) = sum(config.Nsim .* result for (config, result) in zip(configs, results)) ./ mcConfig.Nsim

function average_on_batches(f, mcConfig::BatchedConfig)
    configs = batch_configs(mcConfig)
    return weighted_average(configs, batch_results(f, configs, mcConfig.mode), mcConfig)
end

# The paths of all the batches, in the order of the batches
join_paths(paths::Vector{<:AbstractMatrix}) = reduce(vcat, paths)
join_paths(paths::Vector{<:Vector{<:AbstractMatrix}}) = [reduce(vcat, [batch[i] for batch in paths]) for i in eachindex(paths[1])]

joined_on_batches(f, mcConfig::BatchedConfig) = join_paths(batch_results(f, batch_configs(mcConfig), mcConfig.mode))

simulate(mcProcess::AbstractMonteCarloProcess, und::AbstractUnderlying, zeroCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, T::Number) = joined_on_batches(config -> simulate(mcProcess, und, zeroCurve, config, T), mcConfig)
simulate(mcProcess::VectorialMonteCarloProcess, und::UnderlyingTuple, zeroCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, T::Number) = joined_on_batches(config -> simulate(mcProcess, und, zeroCurve, config, T), mcConfig)
simulate(mcProcess::AbstractMonteCarloEngine, mcConfig::BatchedConfig, T::Number) = joined_on_batches(config -> simulate(mcProcess, config, T), mcConfig)

for (model_type, und_type, payoff_type) in ((BaseProcess, UnderlyingArg, AbstractPayoff), (BaseProcess, UnderlyingArg, Position), (BaseProcess, UnderlyingArg, Array{<:AbstractPayoff}), (BaseProcess, UnderlyingArg, AbstractSmile), (VectorialMonteCarloProcess, UnderlyingTuple, Array{<:Position}))
    @eval function pricer(mcProcess::$model_type, und::$und_type, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, abstractPayoff::$payoff_type)
        return average_on_batches(config -> pricer(mcProcess, und, rfCurve, config, abstractPayoff), mcConfig)
    end
end
function pricer(market::MarketDataSet, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, portfolio::Portfolio)
    return average_on_batches(config -> pricer(market, rfCurve, config, portfolio), mcConfig)
end

# The main run of the control variates is what goes on the batches
function pricer(mcProcess::BaseProcess, und::UnderlyingArg, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:ControlVariates, <:AbstractBatchedMode}, abstractPayoff::AbstractPayoff)
    return pricer_control_variates(mcProcess, und, rfCurve, mcConfig, abstractPayoff)
end

# The means and the variances of the paths of all the batches, without bringing the paths together
function mean_and_variance(mcProcess::BaseProcess, und::UnderlyingArg, rfCurve::AbstractZeroRateCurve, mcConfig::BatchedConfig, abstractPayoffs)
    configs = batch_configs(mcConfig)
    moments = batch_results(config -> payoff_moments(mcProcess, und, rfCurve, config, abstractPayoffs), configs, mcConfig.mode)
    means = weighted_average(configs, first.(moments), mcConfig)
    # each batch holds sample_count(config) independent samples (whole antithetic pairs, see batch_sizes)
    squares = sum(muladd.(sample_count(config), abs2.(moment[1] .- means), moment[2]) for (config, moment) in zip(configs, moments))
    return means, per_path_variance(squares, mcConfig)
end
