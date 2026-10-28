"""
    KouProcess(σ, λ, p, λ₊, λ₋, underlying)

Kou jump diffusion: the log price is a Brownian motion with volatility `σ` plus a compound
Poisson process with intensity `λ` and double exponential jumps (up with probability `p`, with
exponential sizes of rate `λ₊` up and `λ₋` down). The drift is chosen to make the discounted price a martingale.

# Arguments
- `σ`: volatility, must be positive.
- `λ`: jump intensity, must be positive.
- `p`: probability of an up jump, in ``[0, 1]``.
- `λ₊`: rate of the exponential up jumps, must be greater than 1 (for a finite forward).
- `λ₋`: rate of the exponential down jumps, must be positive.
- `underlying`: an [`Underlying`](@ref) with spot price and dividends.
"""
struct KouProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: FiniteActivityProcess{numtype}
    σ::num
    λ::num1
    p::num2
    λ₊::num3
    λ₋::num4
    underlying::abstrUnderlying
    function KouProcess(σ::num, λ::num1, p::num2, λ₊::num3, λ₋::num4, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
        ChainRulesCore.@ignore_derivatives @assert σ > 0 "volatility must be positive"
        ChainRulesCore.@ignore_derivatives @assert λ > 0 "λ, i.e. jump intensity, must be positive"
        ChainRulesCore.@ignore_derivatives @assert 0 <= p <= 1 "p must be a probability"
        ChainRulesCore.@ignore_derivatives @assert λ₊ > 1 "λ₊ must be greater than 1 (finite expected jump)"
        ChainRulesCore.@ignore_derivatives @assert λ₋ > 0 "λ₋ must be positive"
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
