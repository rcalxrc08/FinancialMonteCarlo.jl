"""
    MertonProcess(σ, λ, μ_jump, σ_jump)

Merton jump diffusion: the log price is a Brownian motion with volatility `σ` plus a compound
Poisson process with intensity `λ` and normally distributed jumps. The drift is chosen to make the
discounted price a martingale.

# Arguments
- `σ`: volatility, must be positive.
- `λ`: jump intensity, must be positive.
- `μ_jump`: mean of the jumps of the log price.
- `σ_jump`: standard deviation of the jumps of the log price, must be positive.
"""
struct MertonProcess{num <: Number, num1 <: Number, num2 <: Number, num3 <: Number, numtype <: Number} <: FiniteActivityProcess{numtype}
    σ::num
    λ::num1
    μ_jump::num2
    σ_jump::num3
    function MertonProcess(σ::num, λ::num1, μ_jump::num2, σ_jump::num3) where {num <: Number, num1 <: Number, num2 <: Number, num3 <: Number}
        check_merton(σ, λ, σ_jump)
        zero_typed = zero_of_types(num, num1, num2, num3)
        return new{num, num1, num2, num3, typeof(zero_typed)}(σ, λ, μ_jump, σ_jump)
    end
end

function check_merton(σ, λ, σ_jump)
    check_domain(σ > 0, σ, "volatility must be positive")
    check_domain(λ > 0, λ, "jump intensity must be positive")
    check_domain(σ_jump > 0, σ_jump, "jump volatility must be positive")
end
ChainRulesCore.@non_differentiable check_merton(::Any, ::Any, ::Any)

export MertonProcess;
# The random number of a jump (a normal) and its size from it
jump_raw(::MertonProcess, mcBaseData::AbstractMonteCarloConfiguration) = draw_normal(mcBaseData)
jump_raw_type(::MertonProcess) = Float64
@muladd jump_size(mcProcess::MertonProcess, z) = mcProcess.μ_jump + mcProcess.σ_jump * z
compute_drift(mcProcess::MertonProcess) = exp(mcProcess.μ_jump + mcProcess.σ_jump^2 / 2) - 1
# The jump size at the uniform u, by inversion (AntitheticMC and SobolMC)
@muladd jump_quantile(mcProcess::MertonProcess, u) = mcProcess.μ_jump + mcProcess.σ_jump * norminvcdf(u)
