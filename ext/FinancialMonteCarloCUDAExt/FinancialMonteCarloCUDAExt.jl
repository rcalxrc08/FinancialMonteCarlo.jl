module FinancialMonteCarloCUDAExt

using CUDA
using ChainRulesCore
using Random
using FinancialMonteCarlo: AbstractCudaMode, AbstractMonteCarloConfiguration, AbstractMonteCarloMethod, AbstractSmile, AbstractZeroRateCurve, AmericanPayoff, BaseProcess, BermudanPayoff, FiniteActivityProcess, MonteCarloConfiguration, PathDependentPayoff, VectorialMonteCarloProcess, inner_seed!, integral, maturity, payout_strike, values_at_maturity
import FinancialMonteCarlo: add_jump_matrix!, draw_buffer, draw_normals!, get_array_type, get_matrix_type, payoff, set_seed!, subordinator_host

const CudaConfig = MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractCudaMode}

function set_seed!(mcConfig::CudaConfig)
    if mcConfig.init_rng
        inner_seed!(mcConfig.parallelMode.rng, mcConfig.parallelMode.seed_cpu)
        CUDA.seed!(mcConfig.parallelMode.seed_gpu)
    end
end

get_matrix_type(mcConfig::CudaConfig, ::BaseProcess, price) = CuMatrix{typeof(price)}(undef, mcConfig.Nsim, mcConfig.Nstep + 1);
get_array_type(mcConfig::CudaConfig, ::BaseProcess, price) = CuArray{typeof(price)}(undef, mcConfig.Nsim);
get_matrix_type(::CudaConfig, ::VectorialMonteCarloProcess, price) = Array{CuMatrix{typeof(price)}};

# All the simulate! of the package run here through these methods
draw_buffer(::CudaConfig, numtype, N::Integer) = CuArray{numtype}(undef, N)
draw_normals!(Z, ::CudaConfig) = randn!(Z)
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

# One kernel for all the strikes, instead of one per strike. The result stays on the GPU.
function payoff(S::CuMatrix{num}, smile::AbstractSmile{T_}, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(smile)) where {T_, num <: Number, num2 <: Number}
    T = smile.T
    ST = values_at_maturity(S, T, T1, mcBaseData)
    K = CuArray(reshape(smile.K, 1, :))
    isCall = CuArray(reshape(smile.isCall, 1, :))
    return payout_strike.(ST, K, isCall, T_) .* exp(-integral(rfCurve.r, T))
end

end
