"""
    HestonProcess(σ, σ₀, λ, κ, ρ, θ)

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
"""
struct HestonProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, num5 <: Number, numtype <: Number} <: ItoProcess{numtype}
    σ::num
    σ₀::num1
    λ::num2
    κ::num3
    ρ::num4
    θ::num5
    function HestonProcess(σ::num, σ₀::num1, λ::num2, κ::num3, ρ::num4, θ::num5) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, num5 <: Number}
        check_heston(σ, σ₀, λ, κ, ρ)
        zero_typed = zero_of_types(num, num1, num2, num3, num4, num5)
        return new{num, num1, num2, num3, num4, num5, typeof(zero_typed)}(σ, σ₀, λ, κ, ρ, θ)
    end
end

function check_heston(σ, σ₀, λ, κ, ρ)
    check_domain(σ₀ > 0, σ₀, "initial volatility must be positive")
    check_domain(σ > 0, σ, "volatility of volatility must be positive")
    check_domain(abs(κ + λ) > 1e-14, κ + λ, "unfeasible parameters")
    check_domain(-1 <= ρ <= 1, ρ, "ρ must be a correlation")
end
ChainRulesCore.@non_differentiable check_heston(::Any, ::Any, ::Any, ::Any, ::Any)

export HestonProcess;

# The floor of the variance. Positive, since the derivative of sqrt is infinite at 0, and far from the subnormal numbers:
# eps(0.0) is one, and every operation on a subnormal is about a hundred times slower, which made the paths whose variance
# hits the floor (Feller condition violated) the slowest part of the simulation. sqrt(floatmin) keeps v * dt and sqrt(v)
# normal too, in Float32 and in Float64.
variance_floor(zero_vol) = zero_vol + sqrt(floatmin(eps(zero_vol)))

# The risk neutral parameters, the time step and the zero of the variance, which carries the derivatives of every input
function heston_setup(mcProcess::HestonProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration, T::Number)
    check_maturity(T)
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    κ_s = κ + λ1
    θ_s = κ * mcProcess.θ / κ_s
    dt = T / mcBaseData.Nstep
    drifts = step_drifts(rfCurve.r - dividend(und), dt, mcBaseData.Nstep)
    isDualZeroVol = zero(T * mcProcess.σ₀ * κ_s * θ_s * λ1 * mcProcess.σ * mcProcess.ρ * drift_of_step(drifts, 1))
    return κ_s, θ_s, drifts, dt, isDualZeroVol
end

# One Euler step of the log price and of the variance v (in place) on all the paths of one side, a = ±√dt, with the
# normal of the variance correlated on the fly: e1 ρ + e2_rho √(1 - ρ²). A function of its own, so that the automatic
# differentiation extensions can give it a rule.
function heston_step!(X_to, X_from, v, e1, e2_rho, drift, dt, a, κ, θ, σ, ρ, ρ_orthogonal, v_floor)
    @. X_to = heston_log_step(X_from, v, e1, drift, dt, a)
    @. v = max(heston_variance_step(v, e1, e2_rho, dt, a, κ, θ, σ, ρ, ρ_orthogonal), v_floor)
    return nothing
end
# The two updates on one path, before the floor of the variance (shared with the rules of the AD extensions)
@inline @muladd heston_log_step(x, v, e1, drift, dt, a) = x + drift - v * (dt / 2) + sqrt(v) * a * e1
@inline @muladd heston_variance_step(v, e1, e2_rho, dt, a, κ, θ, σ, ρ, ρ_orthogonal) = v + κ * (θ - v) * dt + σ * sqrt(v) * a * (e1 * ρ + e2_rho * ρ_orthogonal)

function simulate!(X, mcProcess::HestonProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration, T::Number)
    S0 = und.S0
    σ = mcProcess.σ
    ρ = mcProcess.ρ
    κ_s, θ_s, drifts, dt, isDualZeroVol = heston_setup(mcProcess, und, rfCurve, mcBaseData, T)
    isDualZero = isDualZeroVol * S0
    fill!(view(X, :, 1), isDualZero)
    v_floor = variance_floor(isDualZeroVol)
    v_0 = mcProcess.σ₀^2 + isDualZero
    sides = path_sides(X, mcBaseData)
    N = draws_per_side(sides)
    sides_v = map(side -> (side[1], fill!(draw_buffer(mcBaseData, typeof(v_0), N), v_0), side[2]), sides)
    e1 = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero)), N)
    e2_rho = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero)), N)
    ρ_orthogonal = sqrt(1 - ρ * ρ)
    sqrt_dt = sqrt(dt)
    for j = 1:mcBaseData.Nstep
        draw_normals!(e1, mcBaseData)
        draw_normals!(e2_rho, mcBaseData)
        drift = drift_of_step(drifts, j)
        foreach(sides_v) do (X_, v_m, sign_)
            heston_step!(view(X_, :, j + 1), view(X_, :, j), v_m, e1, e2_rho, drift, dt, sign_ * sqrt_dt, κ_s, θ_s, σ, ρ, ρ_orthogonal, v_floor)
        end
    end
    ## Conclude
    exp_paths!(X, S0)
    return
end
