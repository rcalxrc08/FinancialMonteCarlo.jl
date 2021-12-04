"""
Struct for Brownian Motion

		bmProcess=BrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct BrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    σ::num
    μ::CurveType{num1, num4}
    function BrownianMotionVec(σ::num, μ::CurveType{num1, num4}) where {num <: Number, num1 <: Number, num4 <: Number}
        @assert σ > 0 "Volatility must be positive"
        zero_typed = zero(num) + zero(num1) + zero(num4)
        return new{num, num1, num4, typeof(zero_typed)}(σ, μ)
    end
end

function BrownianMotion(σ::num, μ::CurveType{num1, num4}) where {num <: Number, num1 <: Number, num4 <: Number}
    return BrownianMotionVec(σ, μ)
end

function simulate!(X, mcProcess::BrownianMotionVec, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: StandardMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    view(X, :, 1) .= isDualZero
    @inbounds for j = 1:Nstep
        tmp_ = incremental_integral(μ, (j - 1) * dt, dt)
        @inbounds for i = 1:Nsim
            x_i_j = @views X[i, j]
            @views X[i, j+1] = x_i_j + tmp_ + stddev_bm * randn(mcBaseData.rng)
        end
    end

    nothing
end

function simulate!(X, mcProcess::BrownianMotionVec, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    view(X, :, 1) .= isDualZero
    Nsim_2 = div(Nsim, 2)

    @inbounds for j = 1:Nstep
        tmp_ = incremental_integral(μ, (j - 1) * dt, dt)
        @inbounds for i = 1:Nsim_2
            Z = stddev_bm * randn(mcBaseData.rng)
            X[2*i-1, j+1] = X[2*i-1, j] + tmp_ + Z
            X[2*i, j+1] = X[2*i, j] + tmp_ - Z
        end
    end

    nothing
end
