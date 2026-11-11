module FinancialMonteCarloCUDAExt

using CUDA
using ChainRulesCore
using Random
using Adapt
using Distributions: Gamma, InverseGaussian
using FinancialMonteCarlo: AbstractCudaMode, AbstractMonteCarloConfiguration, AbstractSmile, AbstractZeroRateCurve, AmericanPayoff, AsianFixedStrikeOption, AsianFloatingStrikeOption, BarrierPayoff, BermudanPayoff, CudaMode, FiniteActivityProcess, HostAdaptor, JumpList, MonteCarloConfiguration, PathDependentPayoff, add_jumps!, get_rng_type, host_rng, jumps_by_path, maturity, payoff_by_rows, smile_payoffs, to_mode
import FinancialMonteCarlo: add_jump_matrix!, array_type, draw_normals!, draw_subordinator!, mode_adaptor, payoff, seed_device!

# Everything the GPU changes, by the mode: where the arrays live and the generator to seed.
array_type(::AbstractCudaMode) = CuArray
# CuArray and not cu, that would convert to Float32
mode_adaptor(::AbstractCudaMode) = CuArray
Adapt.adapt_storage(::HostAdaptor, x::CuArray) = Array(x)
seed_device!(::AbstractCudaMode, seed::Integer) = CUDA.seed!(seed)

# And by the type of the arrays: the normals of a CuArray come from CURAND.
draw_normals!(Z::CuArray, ::MonteCarloConfiguration) = randn!(Z)

# One thread per element of x, which draws it with the generator of the device. The generator is seeded at every launch
# from the host generator of the configuration, so that set_seed! makes the draws reproducible (CUDA.seed! does not
# reach the generator of the kernels).
function device_sampler_kernel!(x, law, seed::UInt32, ::Type{U}) where {U}
    Random.seed!(Random.default_rng(), seed)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= length(x)
        @inbounds x[i] = device_draw(law, U)
    end
    return
end
function device_draws!(x::CuArray, law, mcBaseData::MonteCarloConfiguration)
    seed = rand(host_rng(mcBaseData), UInt32)
    U = typeof(get_rng_type(zero(eltype(x))))
    threads = 256
    @cuda threads = threads blocks = cld(length(x), threads) device_sampler_kernel!(x, law, seed, U)
    return x
end

# Marsaglia and Tsang (2000): a gamma variate of shape α >= 1 by rejection, with acceptance above 95%; for α < 1 the
# variate of shape α + 1 times U^(1 / α).
function device_draw(law::Gamma, ::Type{U}) where {U}
    α = law.α
    boost = α < 1
    a = boost ? α + 1 : α
    d = a - 1 // 3
    c = 1 / sqrt(9 * d)
    g = zero(d)
    while true
        x = randn(U)
        v = 1 + c * x
        v <= 0 && continue
        v = v * v * v
        u = rand(U)
        x2 = x * x
        if u < 1 - 331 // 10000 * x2 * x2 || log(u) < x2 / 2 + d * (1 - v + log(v))
            g = d * v
            break
        end
    end
    if boost
        g *= rand(U)^(1 / α)
    end
    return g * law.θ
end
# Michael, Schucany and Haas (1976): an inverse Gaussian variate from one normal and one uniform.
function device_draw(law::InverseGaussian, ::Type{U}) where {U}
    μ = law.μ
    λ = law.λ
    ν = randn(U)
    y = ν * ν
    x = μ + μ * μ * y / (2 * λ) - μ / (2 * λ) * sqrt(4 * μ * λ * y + μ * μ * y * y)
    return rand(U) <= μ / (μ + x) ? x : μ * μ / x
end

# The subordinators of variance gamma and normal inverse Gaussian are drawn on the device; any other law on the host,
# with Distributions, and copied.
const DeviceLaw = Union{Gamma, InverseGaussian}
draw_subordinator!(dt_s::CuArray, law, mcBaseData::MonteCarloConfiguration) = copyto!(dt_s, rand(host_rng(mcBaseData), law, length(dt_s)))
draw_subordinator!(dt_s::CuArray, law::DeviceLaw, mcBaseData::MonteCarloConfiguration) = device_draws!(dt_s, law, mcBaseData)

# The jump times and sizes are drawn on the host as a list, sorted by path, and added by one kernel with a thread per
# element of X: it sums the jumps of its path up to its column, so that the writes are contiguous (a thread per path
# writing along its row was 4 times slower).
function jumps_kernel!(X, offsets, steps, sizes)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    Nsim = size(X, 1)
    if idx <= length(X)
        i = (idx - 1) % Nsim + 1
        j = (idx - 1) ÷ Nsim + 1
        @inbounds begin
            first_jump = offsets[i]
            last_jump = offsets[i+1] - 1
            if first_jump <= last_jump
                jumps = zero(eltype(sizes))
                for k = first_jump:last_jump
                    jumps += ifelse(steps[k] <= j, sizes[k], zero(eltype(sizes)))
                end
                X[idx] += jumps
            end
        end
    end
    return
end
function add_jump_matrix!(X::CuMatrix, mcProcess::FiniteActivityProcess, mcConfig::AbstractMonteCarloConfiguration, T::Number)
    # room for the expected number of jumps and 4 standard deviations
    expected = mcConfig.Nsim * mcProcess.λ * T
    jumps = JumpList{eltype(X)}(ceil(Int, expected + 4 * sqrt(expected)))
    add_jumps!(jumps, mcProcess, mcConfig, T)
    isempty(jumps.path) && return
    offsets, steps, sizes = jumps_by_path(jumps, mcConfig.Nsim)
    threads = 256
    @cuda threads = threads blocks = cld(length(X), threads) jumps_kernel!(X, CuArray(offsets), CuArray(steps), CuArray(sizes))
    return
end

# Barrier and Asian payoffs are reductions of the rows of the paths, which run on the device (a new Asian payoff needs its
# asian_payout to join them). The other path dependent
# engines index the paths, and the Longstaff-Schwartz regressions need the host.
for payoff_type in (PathDependentPayoff, Union{AmericanPayoff, BermudanPayoff})
    @eval payoff(S::CuMatrix{num}, payoff_::$payoff_type, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number} = payoff(collect(S), payoff_, rfCurve, mcBaseData, T1)
end
payoff(S::CuMatrix{num}, payoff_::Union{BarrierPayoff, AsianFixedStrikeOption, AsianFloatingStrikeOption}, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number} = payoff_by_rows(S, payoff_, rfCurve, mcBaseData, T1)

# The smile moves to the device first (to_mode), so that a smile with its strikes on the host is priced on the GPU. The
# result stays on the GPU.
payoff(S::CuMatrix, smile::AbstractSmile, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number = maturity(smile)) = smile_payoffs(S, to_mode(smile, CudaMode()), rfCurve, mcBaseData, T1)

end
