"""
    variance(mcProcess, rfCurve, mcConfig, payoff)

Variance of the discounted payoff per path (a vector for a vector of payoffs or a smile), such that the variance of the
price estimate is this value divided by the number of paths `Nsim`. With [`FinancialMonteCarlo.StandardMC`](@ref) it is
the sample variance of the paths; with [`FinancialMonteCarlo.AntitheticMC`](@ref) the paths of a pair `(Z, -Z)` are
correlated, so it is twice the sample variance of the averages of the `Nsim / 2` pairs (#76): smaller than the
variance of the paths for a monotone payoff, larger for an even one (a straddle), for which antithetic sampling
increases the error. With [`FinancialMonteCarlo.SobolMC`](@ref) the points are not independent and this is not an
error estimate.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs, a smile or a [`FinancialMonteCarlo.Position`](@ref).
"""
function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile, Position})
    return mean_and_variance(mcProcess, rfCurve, mcConfig, abstractPayoffs)[2]
end

# The independent samples of the price among the discounted payoffs (one row per path): the paths, or with antithetic
# variates the averages of the pairs, path i and path i + Nsim / 2 (the two sides of path_sides), which are correlated
# (#76). Their number is sample_count.
independent_samples(values, ::MonteCarloConfiguration) = values
function independent_samples(values::AbstractVecOrMat, mcConfig::MonteCarloConfiguration{<:AntitheticMC})
    half = div(mcConfig.Nsim, 2)
    return (selectdim(values, 1, 1:half) .+ selectdim(values, 1, (half+1):(2*half))) ./ 2
end
sample_count(mcConfig::MonteCarloConfiguration) = div(mcConfig.Nsim, paths_per_draw(mcConfig.sampling))

# Mean and sum of the squared deviations of the independent samples, one per payoff: what pools exactly over the
# batches of a parallel mode
path_moments(values, mcConfig::MonteCarloConfiguration) = path_moments(independent_samples(values, mcConfig))
function path_moments(values::AbstractVector)
    mean_ = mean(values)
    return mean_, sum(abs2, values .- mean_)
end

function path_moments(values::AbstractMatrix)
    means = mean(values, dims = 1)
    return vec(means), vec(sum(abs2, values .- means, dims = 1))
end

function payoff_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, payoff_::Union{AbstractPayoff, AbstractSmile})
    abstractPayoff = to_mode(payoff_, execution_mode(mcConfig))
    set_seed!(mcConfig)
    S = simulate(mcProcess, rfCurve, mcConfig, maturity(abstractPayoff))
    return path_moments(payoff(S, abstractPayoff, rfCurve, mcConfig), mcConfig)
end

# A position: the moments of its discounted payoff, the legs combined path by path
function payoff_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, position::Position)
    T_price = price_type(mcProcess, rfCurve, position)
    isempty(position) && return (zero(T_price), zero(T_price))
    set_seed!(mcConfig)
    maxT = maturity(position)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    moments::Tuple{T_price, T_price} = path_moments(payoff(S, position, rfCurve, mcConfig, maxT), mcConfig)
    return moments
end

function payoff_moments(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{<:AbstractPayoff, N}) where {N}
    set_seed!(mcConfig)
    maxT = maturity(abstractPayoffs)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    # An array of payoffs of different types has an abstract element type: the type of the moments comes from
    # price_type, as in pricer, so that the result is inferred
    T_price = price_type(mcProcess, rfCurve, abstractPayoffs)
    moments::Array{Tuple{T_price, T_price}, N} = [path_moments(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT), mcConfig) for abstractPayoff in abstractPayoffs]
    return first.(moments), last.(moments)
end

function mean_and_variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs)
    means, squares = payoff_moments(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    return means, per_path_variance(squares, mcConfig)
end

# The sum of the squared deviations of the n independent samples to the variance per path: the sample variance of the
# samples times the paths per sample (integers: Float32 stays Float32)
per_path_variance(squares, mcConfig::MonteCarloConfiguration) = squares .* paths_per_draw(mcConfig.sampling) ./ (sample_count(mcConfig) - 1)
