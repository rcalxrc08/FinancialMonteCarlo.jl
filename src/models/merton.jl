"""
Struct for Merton Process

		mertonProcess=MertonProcess(σ::num1,λ::num2,μ_jump::num3,σ_jump::num4) where {num1,num2,num3,num4<: Number}
	
Where:

		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		μ_jump =	jumps mean.
		σ_jump =	jumps variance.
"""
struct MertonProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: FiniteActivityProcess{numtype}
    σ::num
    λ::num1
    μ_jump::num2
    σ_jump::num3
    underlying::abstrUnderlying
    function MertonProcess(σ::num, λ::num1, μ_jump::num2, σ_jump::num3, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying}
        ChainRulesCore.@ignore_derivatives @assert σ > 0 "volatility must be positive"
        ChainRulesCore.@ignore_derivatives @assert λ > 0 "jump intensity must be positive"
        ChainRulesCore.@ignore_derivatives @assert σ_jump > 0 "positive λ must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num2) + zero(num3)
        return new{num, num1, num2, num3, abstrUnderlying, typeof(zero_typed)}(σ, λ, μ_jump, σ_jump, underlying)
    end
end

export MertonProcess;
compute_jump_size(mcProcess::MertonProcess, mcBaseData::MonteCarloConfiguration) = mcProcess.μ_jump + mcProcess.σ_jump * randn(mcBaseData.parallelMode.rng);
compute_drift(mcProcess::MertonProcess) = exp(mcProcess.μ_jump + mcProcess.σ_jump^2 / 2.0) - 1.0
