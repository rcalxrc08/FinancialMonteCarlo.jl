"""
    VarianceGammaProcess(σ, θ, κ)

Variance Gamma model: the log price is a Brownian motion with drift `θ` and volatility `σ`,
evaluated at a gamma subordinator with unit mean rate and variance rate `κ`. The drift is chosen to
make the discounted price a martingale.

# Arguments
- `σ`: volatility of the subordinated Brownian motion, must be positive.
- `θ`: drift of the subordinated Brownian motion.
- `κ`: variance rate of the subordinator, must be positive.

The parameters must satisfy ``1 - σ^2 κ / 2 - θ κ > 0`` for the forward to be finite.
"""
struct VarianceGammaProcess{num <: Number, num1 <: Number, num2 <: Number, numtype <: Number} <: InfiniteActivityProcess{numtype}
    σ::num
    θ::num1
    κ::num2
    function VarianceGammaProcess(σ::num, θ::num1, κ::num2) where {num <: Number, num1 <: Number, num2 <: Number}
        check_vg(σ, θ, κ)
        zero_typed = zero_of_types(num, num1, num2)
        return new{num, num1, num2, typeof(zero_typed)}(σ, θ, κ)
    end
end

function check_vg(σ, θ, κ)
    check_domain(σ > 0, σ, "volatility must be positive")
    check_domain(κ > 0, κ, "κappa must be positive")
    check_domain(1 - σ * σ * κ / 2 - θ * κ >= 0, (σ, θ, κ), "Parameters with unfeasible values")
end
ChainRulesCore.@non_differentiable check_vg(::Any, ::Any, ::Any)

export VarianceGammaProcess;

@muladd compute_drift(::VarianceGammaProcess, σ::Number, θ::Number, κ::Number) = -1 / κ * log(1 - σ^2 * κ / 2 - θ * κ)
compute_drift(mcProcess::VarianceGammaProcess) = compute_drift(mcProcess, mcProcess.σ, mcProcess.θ, mcProcess.κ)
compute_subordinator_law(::VarianceGammaProcess, κ::Number, dt::Number) = Gamma(dt / κ, κ)
compute_subordinator_law(mcProcess::VarianceGammaProcess, dt::Number) = compute_subordinator_law(mcProcess, mcProcess.κ, dt)
