"""
    distribution(mcProcess, underlying, rfCurve, mcConfig, T)

Simulated values of the underlying at time `T`, one per path (the last column of the paths
returned by [`simulate`](@ref)).

# Arguments
- `mcProcess`: the model.
- `underlying`: the [`Underlying`](@ref) (spot and dividends), a tuple of them for a multivariate model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `T`: time, must be positive.
"""
function distribution(mcProcess::BaseProcess, und::UnderlyingArg, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, T::num_) where {num_ <: Number}
    set_seed!(mcConfig)
    S = simulate(mcProcess, und, rfCurve, mcConfig, T)

    return S[:, end]
end
