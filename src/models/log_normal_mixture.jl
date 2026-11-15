"""
    LogNormalMixture(η, λ, underlying)

LogNormal mixture model of Brigo and Mercurio: the local volatility diffusion, in ``X_t = \\log(S_t / S_0)``,
```math
dX_t = \\left(r_t - d_t - \\frac{σ^2(t, X_t)}{2}\\right) dt + σ(t, X_t)\\,dW_t, \\qquad
σ^2(t, x) = \\frac{\\sum_i λ_i η_i^2 p_i(t, x)}{\\sum_i λ_i p_i(t, x)},
```
where ``p_i(t, \\cdot)`` is the density of ``X_t`` in the Black-Scholes model with volatility ``η_i``. The law of
``S_t`` is the mixture ``\\sum_i λ_i \\text{Law}(S^i_t)`` at every time, so European prices are the weighted sums of the
Black-Scholes prices with volatilities ``η_i``.

Simulated with an Euler scheme on ``X``, the local volatility taken at the start of each step (at ``t = 0``, its limit
``\\sum_i λ_i η_i / \\sum_i λ_i / η_i``): European prices converge to the closed form as the number of steps grows,
slowly when the volatilities are far apart (the local volatility then changes sharply with ``X``).

# Examples
```julia
LogNormalMixture((0.2, 0.3), (0.4,), Underlying(100.0, 0.01))   # weights 0.4 and 0.6
```

# Arguments
- `η`: tuple of the ``n`` volatilities, all positive.
- `λ`: tuple of the first ``n - 1`` weights, positive and summing to at most 1 (the last weight is `1 - sum(λ)`).
- `underlying`: an [`Underlying`](@ref) with spot price and dividends.

`η` and `λ` are tuples, so the number of components is part of the type: a model can be passed to a GPU kernel, and
the loops over the components are unrolled.
"""
struct LogNormalMixture{N, M, num <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ItoProcess{numtype}
    η::NTuple{N, num}
    λ::NTuple{M, num2}
    underlying::abstrUnderlying
    # at least one weight, so two components: the number types are those of the first elements
    function LogNormalMixture(η::Tuple{num, Vararg{num}}, λ::Tuple{num2, Vararg{num2}}, underlying::abstrUnderlying) where {num <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying}
        N, M = length(η), length(λ)
        ChainRulesCore.@ignore_derivatives @assert M == N - 1 "length(λ) must be length(η) - 1"
        ChainRulesCore.@ignore_derivatives @assert minimum(η) > 0 "Volatilities must be positive"
        ChainRulesCore.@ignore_derivatives @assert minimum(λ) > 0 "weights must be positive"
        ChainRulesCore.@ignore_derivatives @assert sum(λ) <= 1 "λs must be weights"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num2)
        return new{N, M, num, num2, abstrUnderlying, typeof(zero_typed)}(η, λ, underlying)
    end
end

export LogNormalMixture;

# All the weights, the last one 1 - sum(λ)
mixture_weights(model::LogNormalMixture) = (model.λ..., 1 - sum(model.λ))

# The local variance at x = log(S / F_t), F_t the forward, t > 0. The weight of component i is λ_i p_i(t, x), and the
# densities share their factor 1 / √(2πt), so it is exp(e_i): the log of λ_i / η_i minus the square of the distance of
# x from the mean of the component, -η_i² t / 2, over twice its variance η_i² t. The exponents are shifted by their
# maximum (log-sum-exp), so that in the far tails at small t the weights do not all underflow to 0 / 0.
function lnm_local_variance(x, t, η::NTuple, logw::NTuple)
    e = map((η_i, logw_i) -> logw_i - (x + η_i^2 * t / 2)^2 / (2 * η_i^2 * t), η, logw)
    e_max = maximum(e)
    w = map(e_i -> exp(e_i - e_max), e)
    return sum(map((w_i, η_i) -> w_i * η_i^2, w, η)) / sum(w)
end

# At t = 0 every path is at the spot, and the weights tend to λ_i / η_i: σ²(0, 0) = Σ λ_i η_i / Σ λ_i / η_i
lnm_initial_variance(η::NTuple, λ::NTuple) = sum(map(*, λ, η)) / sum(map(/, λ, η))

function lnm_setup(mcProcess::LogNormalMixture, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0
    η = mcProcess.η
    λ = mixture_weights(mcProcess)
    logw = map((λ_i, η_i) -> log(λ_i) - log(η_i), λ, η)
    μ = rfCurve.r - dividend(mcProcess)
    return η, logw, lnm_initial_variance(η, λ), μ, T / mcBaseData.Nstep
end

# One Euler substep of the log price x from time t, of length h, with the drift of the log price over it and the normal z
@inline function lnm_euler(x, t, h, drift, z, carry, σ2_0, η, logw)
    σ2 = iszero(t) ? σ2_0 : lnm_local_variance(x - carry, t, η, logw)
    return x + drift - σ2 * h / 2 + sqrt(σ2 * h) * z
end

function simulate!(X, mcProcess::LogNormalMixture, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration, T::Number)
    η, logw, σ2_0, μ, dt = lnm_setup(mcProcess, rfCurve, mcBaseData, T)
    S0 = mcProcess.underlying.S0
    fill!(view(X, :, 1), zero(eltype(X)))
    sides = path_sides(X, mcBaseData)
    Z = draw_buffer(mcBaseData, typeof(get_rng_type(zero(eltype(X)))), draws_per_side(sides))
    drifts = step_drifts(μ, dt, mcBaseData.Nstep)
    for j = 1:mcBaseData.Nstep
        draw_normals!(Z, mcBaseData)
        t = (j - 1) * dt
        drift = drift_of_step(drifts, j)
        carry = integral(μ, t)
        foreach(sides) do (X_, sign_)
            @views @. X_[:, j+1] = lnm_euler(X_[:, j], t, dt, drift, sign_ * Z, carry, σ2_0, (η,), (logw,))
        end
    end
    @. X = S0 * exp(X)
    return nothing
end

function simulate!(X, mcProcess::LogNormalMixture, rfCurve::AbstractZeroRateCurve, mcBaseData::SobolConfig, T::Number)
    η, logw, σ2_0, μ, dt = lnm_setup(mcProcess, rfCurve, mcBaseData, T)
    S0 = mcProcess.underlying.S0
    Nstep = mcBaseData.Nstep
    fill!(view(X, :, 1), zero(eltype(X)))
    drifts = step_drifts(μ, dt, Nstep)
    carries = [integral(μ, (j - 1) * dt) for j = 1:Nstep]
    on_sobol_paths(Nstep, mcBaseData, typeof(get_rng_type(zero(eltype(X))))) do i, normals
        @inbounds for j = 1:Nstep
            X[i, j+1] = lnm_euler(X[i, j], (j - 1) * dt, dt, drift_of_step(drifts, j), normals[j], carries[j], σ2_0, η, logw)
        end
    end
    @. X = S0 * exp(X)
    return nothing
end
