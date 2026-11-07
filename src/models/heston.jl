"""
    HestonProcess(σ, σ₀, λ, κ, ρ, θ, underlying)

Heston stochastic volatility model:
```math
dS_t = (r_t - d_t) S_t\\,dt + \\sqrt{v_t} S_t\\,dW_t, \\qquad
dv_t = κ^*(θ^* - v_t)\\,dt + σ \\sqrt{v_t}\\,dZ_t, \\qquad d\\langle W, Z\\rangle_t = ρ\\,dt,
```
with risk neutral parameters ``κ^* = κ + λ`` and ``θ^* = κθ / (κ + λ)``. Both are simulated with an Euler
scheme on the log price, flooring the variance at a tiny positive value.

# Arguments
- `σ`: volatility of the variance, must be positive.
- `σ₀`: initial volatility, must be positive.
- `λ`: market price of volatility risk.
- `κ`: speed of mean reversion (``κ + λ`` must not be zero).
- `ρ`: correlation between the Brownian motions, in ``[-1, 1]``.
- `θ`: long term variance.
- `underlying`: an [`Underlying`](@ref) with spot price and dividends.
"""
struct HestonProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, num5 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ItoProcess{numtype}
    σ::num
    σ₀::num1
    λ::num2
    κ::num3
    ρ::num4
    θ::num5
    underlying::abstrUnderlying
    function HestonProcess(σ::num, σ₀::num1, λ::num2, κ::num3, ρ::num4, θ::num5, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, num5 <: Number, abstrUnderlying <: AbstractUnderlying}
        ChainRulesCore.@ignore_derivatives @assert σ₀ > 0.0 "initial volatility must be positive"
        ChainRulesCore.@ignore_derivatives @assert σ > 0.0 "volatility of volatility must be positive"
        ChainRulesCore.@ignore_derivatives @assert abs(κ + λ) > 1e-14 "unfeasible parameters"
        ChainRulesCore.@ignore_derivatives @assert (-1.0 <= ρ <= 1.0) "ρ must be a correlation"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num2) + zero(num3) + zero(num4) + zero(num5)
        return new{num, num1, num2, num3, num4, num5, abstrUnderlying, typeof(zero_typed)}(σ, σ₀, λ, κ, ρ, θ, underlying)
    end
end

export HestonProcess;

# The floor of the variance. Positive, since the derivative of sqrt is infinite at 0, and far from the subnormal numbers:
# eps(0.0) is one, and every operation on a subnormal is about a hundred times slower, which made the paths whose variance
# hits the floor (Feller condition violated) the slowest part of the simulation. sqrt(floatmin) keeps v * dt and sqrt(v)
# normal too, in Float32 and in Float64.
variance_floor(zero_vol) = zero_vol + sqrt(floatmin(eps(zero_vol)))

# The risk neutral parameters, the time step and the zero of the variance, which carries the derivatives of every input
function heston_setup(mcProcess::HestonProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    κ_s = κ + λ1
    θ_s = κ * mcProcess.θ / κ_s
    dt = T / mcBaseData.Nstep
    drifts = step_drifts(rfCurve.r - dividend(mcProcess), dt, mcBaseData.Nstep)
    isDualZeroVol = zero(T * mcProcess.σ₀ * κ_s * θ_s * λ1 * mcProcess.σ * mcProcess.ρ * drift_of_step(drifts, 1))
    return κ_s, θ_s, drifts, dt, isDualZeroVol
end

function simulate!(X, mcProcess::HestonProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration, T::Number)
    S0 = mcProcess.underlying.S0
    σ = mcProcess.σ
    ρ = mcProcess.ρ
    κ_s, θ_s, drifts, dt, isDualZeroVol = heston_setup(mcProcess, rfCurve, mcBaseData, T)
    isDualZero = isDualZeroVol * S0
    fill!(view(X, :, 1), isDualZero)
    v_floor = variance_floor(isDualZeroVol)
    v_0 = mcProcess.σ₀^2 + isDualZero
    sides = path_sides(X, mcBaseData)
    N = draws_per_side(sides)
    sides_v = map(side -> (side[1], fill!(draw_buffer(mcBaseData, typeof(v_0), N), v_0), side[2]), sides)
    e1 = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero)), N)
    e2_rho = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero)), N)
    e2 = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero) + zero(ρ)), N)
    ρ_orthogonal = sqrt(1 - ρ * ρ)
    sqrt_dt = sqrt(dt)
    for j = 1:mcBaseData.Nstep
        draw_normals!(e1, mcBaseData)
        draw_normals!(e2_rho, mcBaseData)
        @. e2 = e1 * ρ + e2_rho * ρ_orthogonal
        drift = drift_of_step(drifts, j)
        for (X_, v_m, sign_) in sides_v
            sqrt_dt_step = sign_ * sqrt_dt
            @views @. X_[:, j+1] = X_[:, j] + (drift - v_m * dt / 2) + sqrt(v_m) * sqrt_dt_step * e1
            @. v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt_dt_step * e2
            @. v_m = max(v_m, v_floor)
        end
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end
