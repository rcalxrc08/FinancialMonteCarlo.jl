"""
    variance(mcProcess, rfCurve, mcConfig, payoff)

Sample variance of the discounted payoff over the simulated paths (a vector for a vector of payoffs or a smile).
The variance of the price estimate is this value divided by the number of paths.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs or a smile.
"""
function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile})
    return mean_and_variance(mcProcess, rfCurve, mcConfig, abstractPayoffs)[2]
end

# Mean and sum of the squared deviations of the discounted payoffs, one per payoff: what pools exactly over the batches of a parallel mode
function path_moments(Payoff::AbstractVector)
    mean_ = mean(Payoff)
    return mean_, sum(abs2, Payoff .- mean_)
end

function path_moments(Payoffs::AbstractMatrix)
    means = mean(Payoffs, dims = 1)
    return vec(means), vec(sum(abs2, Payoffs .- means, dims = 1))
end

function payoff_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, payoff_::Union{AbstractPayoff, AbstractSmile})
    abstractPayoff = to_mode(payoff_, execution_mode(mcConfig))
    set_seed!(mcConfig)
    S = simulate(mcProcess, rfCurve, mcConfig, maturity(abstractPayoff))
    return path_moments(payoff(S, abstractPayoff, rfCurve, mcConfig))
end

function payoff_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{<:AbstractPayoff})
    set_seed!(mcConfig)
    maxT = maturity(abstractPayoffs)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    moments = [path_moments(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]
    return first.(moments), last.(moments)
end

function mean_and_variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs)
    Means, squares = payoff_moments(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    return Means, squares ./ (mcConfig.Nsim - 1)
end
