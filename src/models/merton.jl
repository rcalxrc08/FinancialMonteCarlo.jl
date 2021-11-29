"""
Struct for Merton Process

		mertonProcess=MertonProcess(σ::num1,λ::num2,μJ::num3,σJ::num4) where {num1,num2,num3,num4<: Number}
	
Where:\n
		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		μJ =	jumps mean.
		σJ =	jumps variance.
"""
mutable struct MertonProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: FiniteActivityProcess{numtype}
    σ::num
    λ::num1
    μJ::num2
    σJ::num3
    underlying::abstrUnderlying
    function MertonProcess(σ::num, λ::num1, μJ::num2, σJ::num3, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying}
        @assert σ > 0 "volatility must be positive"
        @assert λ > 0 "jump intensity must be positive"
        @assert σJ > 0 "positive λ must be positive"
        zero_typed = zero(num) + zero(num1) + zero(num2) + zero(num3)
        return new{num, num1, num2, num3, abstrUnderlying, typeof(zero_typed)}(σ, λ, μJ, σJ, underlying)
    end
end

export MertonProcess;
compute_jump_size(mcProcess::MertonProcess, mcBaseData::MonteCarloConfiguration) = mcProcess.μJ + mcProcess.σJ * randn(mcBaseData.rng);
compute_drift(mcProcess::MertonProcess) = exp(mcProcess.μJ + mcProcess.σJ^2 / 2.0) - 1.0
