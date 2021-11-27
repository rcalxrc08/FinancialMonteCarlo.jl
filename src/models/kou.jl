"""
Struct for Kou Process

		kouProcess=KouProcess(σ::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		p  =	prob. of having positive jump.
		λ₊ =	positive jump size.
		λ₋ =	negative jump size.
"""
mutable struct KouProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: FiniteActivityProcess{numtype}
    σ::num
    λ::num1
    p::num2
    λ₊::num3
    λ₋::num4
    underlying::abstrUnderlying
    function KouProcess(σ::num, λ::num1, p::num2, λ₊::num3, λ₋::num4, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ <= 0
            error("volatility must be positive")
        elseif λ <= 0
            error("jump intensity must be positive")
        elseif λ₊ <= 0
            error("positive λ must be positive")
        elseif λ₋ <= 0
            error("negative λ must be positive")
        elseif !(0 <= p <= 1)
            error("p must be a probability")
        else
            zero_typed = promote_type(num, num1, num2, num3, num4)
            return new{num, num1, num2, num3, num4, abstrUnderlying, zero_typed}(σ, λ, p, λ₊, λ₋, underlying)
        end
    end
end

export KouProcess;

function compute_jump_size(mcProcess::KouProcess, mcBaseData::MonteCarloConfiguration)
    return if rand(mcBaseData.rng) < mcProcess.p
        quantile_exp(mcProcess.λ₊, rand(mcBaseData.rng))
    else
        -quantile_exp(mcProcess.λ₋, rand(mcBaseData.rng))
    end
end
compute_drift(mcProcess::KouProcess) = mcProcess.p / (mcProcess.λ₊ - 1) - (1 - mcProcess.p) / (mcProcess.λ₋ + 1)
