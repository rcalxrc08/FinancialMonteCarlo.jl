"""
Struct for Kou Process

		kouProcess=KouProcess(σ::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:

		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		p  =	prob. of having positive jump.
		λ₊ =	positive jump size.
		λ₋ =	negative jump size.
"""
struct KouProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: FiniteActivityProcess{numtype}
    σ::num
    λ::num1
    p::num2
    λ₊::num3
    λ₋::num4
    underlying::abstrUnderlying
    function KouProcess(σ::num, λ::num1, p::num2, λ₊::num3, λ₋::num4, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
        @assert σ > 0 "volatility must be positive"
        @assert λ > 0 "λ, i.e. jump intensity, must be positive"
        @assert 0 <= p <= 1 "p must be a probability"
        @assert λ₊ > 0 "λ₊ must be positive"
        @assert λ₋ > 0 "λ₋ must be positive"
        zero_typed = promote_type(num, num1, num2, num3, num4)
        return new{num, num1, num2, num3, num4, abstrUnderlying, zero_typed}(σ, λ, p, λ₊, λ₋, underlying)
    end
end

export KouProcess;

function compute_jump_size(mcProcess::KouProcess, mcBaseData::MonteCarloConfiguration)
    if rand(mcBaseData.parallelMode.rng) < mcProcess.p
        return randexp(mcBaseData.parallelMode.rng) / mcProcess.λ₊
    end
    return -randexp(mcBaseData.parallelMode.rng) / mcProcess.λ₋
end
compute_drift(mcProcess::KouProcess) = mcProcess.p / (mcProcess.λ₊ - 1) - (1 - mcProcess.p) / (mcProcess.λ₋ + 1)
