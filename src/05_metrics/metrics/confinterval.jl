using Distributions

"""
    confinter(mcProcess, rfCurve, mcConfig, payoff, alpha = 99 // 100)

Confidence interval at level `alpha` for the price of `payoff`, based on the Student t approximation
of the Monte Carlo estimator, from the independent samples: the paths, or the antithetic pairs with
[`FinancialMonteCarlo.AntitheticMC`](@ref) (see [`variance`](@ref)). Not an error estimate with
[`FinancialMonteCarlo.SobolMC`](@ref), whose points are not independent. Returns a tuple `(lower, upper)`, or a vector of them for a vector of payoffs or a smile.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs, a smile or a [`FinancialMonteCarlo.Position`](@ref).
- `alpha`: confidence level.
"""
function confinter(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile, Position}, alpha::Real = 99 // 100)
    means, vars = mean_and_variance(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    return student_intervals(means, vars, mcConfig.Nsim, alpha, sample_count(mcConfig))
end

# the mean of n independent samples: n - 1 degrees of freedom
student_quantile(n, alpha) = quantile(Distributions.TDist(n - 1), (1 + alpha) / 2)

# The quantile comes in the precision of its library: the interval stays in the one of the prices
function student_interval(mean1, var1, Nsim, tstar)
    standard_error = sqrt(var1 / Nsim)
    tstar_ = oftype(standard_error, tstar)
    return muladd(-tstar_, standard_error, mean1), muladd(tstar_, standard_error, mean1)
end

# vars per path (see variance): the standard error of the price is sqrt(var / Nsim); n independent samples
student_intervals(mean1::Number, var1::Number, Nsim, alpha, n) = student_interval(mean1, var1, Nsim, student_quantile(n, alpha))

function student_intervals(means, vars, Nsim, alpha, n)
    tstar = student_quantile(n, alpha)
    return [student_interval(mean1, var1, Nsim, tstar) for (mean1, var1) in zip(convert(Array, means), convert(Array, vars))]
end
