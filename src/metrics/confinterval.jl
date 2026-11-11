using Distributions

"""
    confinter(mcProcess, rfCurve, mcConfig, payoff, alpha = 99 // 100)

Confidence interval at level `alpha` for the price of `payoff`, based on the Student t approximation
of the Monte Carlo estimator. Returns a tuple `(lower, upper)`, or a vector of them for a vector of payoffs or a smile.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs, a smile or a [`FinancialMonteCarlo.Position`](@ref).
- `alpha`: confidence level.
"""
function confinter(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile, Position}, alpha::Real = 99 // 100)
    means, vars = mean_and_variance(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    return student_intervals(means, vars, mcConfig.Nsim, alpha)
end

student_quantile(Nsim, alpha) = quantile(Distributions.TDist(Nsim), (1 + alpha) / 2)

# The quantile comes in the precision of its library: the interval stays in the one of the prices
function student_interval(mean1, var1, Nsim, tstar)
    standard_error = sqrt(var1 / Nsim)
    tstar_ = oftype(standard_error, tstar)
    return muladd(-tstar_, standard_error, mean1), muladd(tstar_, standard_error, mean1)
end

student_intervals(mean1::Number, var1::Number, Nsim, alpha) = student_interval(mean1, var1, Nsim, student_quantile(Nsim, alpha))

function student_intervals(means, vars, Nsim, alpha)
    tstar = student_quantile(Nsim, alpha)
    return [student_interval(mean1, var1, Nsim, tstar) for (mean1, var1) in zip(convert(Array, means), convert(Array, vars))]
end
