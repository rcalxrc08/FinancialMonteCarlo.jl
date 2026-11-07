module FinancialMonteCarloCUDAExt

using CUDA
using ChainRulesCore
using Random
using Adapt
using FinancialMonteCarlo: AbstractCudaMode, AbstractMonteCarloConfiguration, AbstractSmile, AbstractZeroRateCurve, AmericanPayoff, BermudanPayoff, CudaMode, FiniteActivityProcess, HostAdaptor, MonteCarloConfiguration, PathDependentPayoff, integral, maturity, payout_strike, to_mode, values_at_maturity
import FinancialMonteCarlo: add_jump_matrix!, array_type, draw_normals!, mode_adaptor, payoff, seed_device!, subordinator_host

# Everything the GPU changes, by the mode: where the arrays live and the generator to seed.
array_type(::AbstractCudaMode) = CuArray
# CuArray and not cu, that would convert to Float32
mode_adaptor(::AbstractCudaMode) = CuArray
Adapt.adapt_storage(::HostAdaptor, x::CuArray) = Array(x)
seed_device!(::AbstractCudaMode, seed::Integer) = CUDA.seed!(seed)

# And by the type of the arrays: the normals of a CuArray come from CURAND, the laws of Distributions sample on the host.
draw_normals!(Z::CuArray, ::MonteCarloConfiguration) = randn!(Z)
subordinator_host(dt_s::CuArray) = Array{eltype(dt_s)}(undef, length(dt_s))

# The jump times are drawn on the host, and writing them to X would launch a kernel for each jump.
function add_jump_matrix!(X::CuMatrix, mcProcess::FiniteActivityProcess, mcConfig::AbstractMonteCarloConfiguration, T::Number)
    X_jumps = zeros(eltype(X), mcConfig.Nsim, mcConfig.Nstep + 1)
    add_jump_matrix!(X_jumps, mcProcess, mcConfig, T)
    X .+= CuArray(X_jumps)
    return
end

# The path dependent engines index the paths, and the Longstaff-Schwartz regressions need the host.
for payoff_type in (PathDependentPayoff, Union{AmericanPayoff, BermudanPayoff})
    @eval payoff(S::CuMatrix{num}, payoff_::$payoff_type, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number} = payoff(collect(S), payoff_, rfCurve, mcBaseData, T1)
end

# One kernel for all the strikes, instead of one per strike. The result stays on the GPU. The metrics move the smile to
# the device once (to_mode), so that here it is already there.
function payoff(S::CuMatrix{num}, smile_::AbstractSmile{T_}, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(smile_)) where {T_, num <: Number, num2 <: Number}
    smile = to_mode(smile_, CudaMode())
    T = smile.T
    ST = values_at_maturity(S, T, T1, mcBaseData)
    K = reshape(smile.K, 1, :)
    isCall = reshape(smile.isCall, 1, :)
    return payout_strike.(ST, K, isCall, T_) .* exp(-integral(rfCurve.r, T))
end

end
