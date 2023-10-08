"""
Struct for Heston Process

		hstProcess=HestonProcess(σ::num1,σ₀::num2,λ::num3,κ::num4,ρ::num5,θ::num6) where {num1,num2,num3,num4,num5,num6 <: Number}
	
Where:

		σ	=	volatility of volatility of the process.
		σ₀	=   initial volatility of the process.
		λ	=	shift value of the process.
		κ	=	mean reversion of the process.
		ρ	=	correlation between brownian motions.
		θ	=	long value of square of volatility of the process.
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

function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRate, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    ChainRulesCore.@ignore_derivatives @assert T > 0.0

    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    isDualZeroVol = zero(T * r * σ₀ * θ_s * κ_s * σ * ρ)
    isDualZero = isDualZeroVol * S0
    view(X, :, 1) .= isDualZero
    v_m = [σ₀^2 + isDualZero for _ = 1:Nsim]
    isDualZero_eps = isDualZeroVol + eps(isDualZeroVol)
    e1 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    #Wrong, is dual zero shouldn't have rho
    e2_rho = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    e2 = Array{typeof(get_rng_type(isDualZero) + zero(ρ))}(undef, Nsim)
    tmp_cost = sqrt(1 - ρ * ρ)
    #TODO: acnaoicna
    for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, e1)
        randn!(mcBaseData.parallelMode.rng, e2_rho)
        @. e2 = e1 * ρ + e2_rho * tmp_cost
        @views @. X[:, j+1] = X[:, j] + ((r - d) - 0.5 * v_m) * dt + sqrt(v_m) * sqrt(dt) * e1
        @. v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt(dt) * e2
        #when v_m = 0.0, the derivative becomes NaN
        @. v_m = max(v_m, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end

function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRate, mcBaseData::SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    ChainRulesCore.@ignore_derivatives @assert T > 0.0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    isDualZeroVol = zero(T * r * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = isDualZeroVol * S0
    view(X, :, 1) .= isDualZero
    isDualZero_eps = isDualZeroVol + eps(isDualZeroVol)
    Nsim_2 = div(Nsim, 2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    v_m_1 = [σ₀^2 + isDualZero for _ = 1:Nsim_2]
    v_m_2 = [σ₀^2 + isDualZero for _ = 1:Nsim_2]
    e1 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2_rho = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, e1)
        randn!(mcBaseData.parallelMode.rng, e2_rho)
        @. e2 = -(e1 * ρ + e2_rho * sqrt(1 - ρ * ρ))
        @views @. Xp[:, j+1] = Xp[:, j] + ((r - d) - 0.5 * v_m_1) * dt + sqrt(v_m_1) * sqrt(dt) * e1
        @views @. Xm[:, j+1] = Xm[:, j] + ((r - d) - 0.5 * v_m_2) * dt - sqrt(v_m_2) * sqrt(dt) * e1
        @. v_m_1 += κ_s * (θ_s - v_m_1) * dt + σ * sqrt(v_m_1) * sqrt(dt) * e2
        @. v_m_2 += κ_s * (θ_s - v_m_2) * dt - σ * sqrt(v_m_2) * sqrt(dt) * e2
        @. v_m_1 = max(v_m_1, isDualZero_eps)
        @. v_m_2 = max(v_m_2, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end
