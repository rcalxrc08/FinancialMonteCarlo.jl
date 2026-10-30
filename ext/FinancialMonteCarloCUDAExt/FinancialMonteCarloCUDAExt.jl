module FinancialMonteCarloCUDAExt

using CUDA
using ChainRulesCore
using Random
using Statistics: quantile
using FinancialMonteCarlo: AbstractCudaMode, AbstractMonteCarloConfiguration, AbstractMonteCarloMethod, AbstractZeroRateCurve, AmericanPayoff, BaseProcess, BermudanPayoff, BrownianMotion, BrownianMotionVec, CudaAntitheticMonteCarloConfig, CudaMode, CudaMonteCarloConfig, FiniteActivityProcess, HestonProcess, MonteCarloConfiguration, PathDependentPayoff, SubordinatedBrownianMotion, SubordinatedBrownianMotionVec, VectorialMonteCarloProcess, ZeroRate, ZeroRateCurve, add_jump_matrix!, compute_drift, dividend, get_rng_type, incremental_integral, inner_seed!, maturity
import FinancialMonteCarlo: get_array_type, get_matrix_type, payoff, set_seed!, simulate!

include("brownian_motion_cuda.jl")
include("brownian_motion_cuda_aug.jl")
include("heston_cuda.jl")
include("heston_aug_cuda.jl")
include("fa_levy_cuda.jl")
include("subordinated_brownian_motion_cuda.jl")
include("subordinated_brownian_motion_aug_cuda.jl")

function set_seed!(mcConfig::MonteCarloConfiguration{type1, type2, type3, type4}) where {type1 <: Integer, type2 <: Integer, type3 <: AbstractMonteCarloMethod, type4 <: AbstractCudaMode}
    if mcConfig.init_rng
        inner_seed!(mcConfig.parallelMode.rng, mcConfig.parallelMode.seed_cpu)
        CUDA.seed!(mcConfig.parallelMode.seed_gpu)
    end
end

get_matrix_type(mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractCudaMode}, ::BaseProcess, price) = CuMatrix{typeof(price)}(undef, mcConfig.Nsim, mcConfig.Nstep + 1);
get_array_type(mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractCudaMode}, ::BaseProcess, price) = CuArray{typeof(price)}(undef, mcConfig.Nsim);
get_matrix_type(::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractCudaMode}, ::VectorialMonteCarloProcess, price) = Array{CuMatrix{typeof(price)}};

function payoff(S::CuMatrix{num}, payoff_::PathDependentPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number}
    S_ = collect(S)

    return payoff(S_, payoff_, rfCurve, mcBaseData, T1)
end

# The Longstaff-Schwartz regressions index the paths, so the exercise engines run on the host.
payoff(S::CuMatrix{num}, payoff_::AmericanPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number} = payoff(collect(S), payoff_, rfCurve, mcBaseData, T1)
payoff(S::CuMatrix{num}, payoff_::BermudanPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number} = payoff(collect(S), payoff_, rfCurve, mcBaseData, T1)

end
