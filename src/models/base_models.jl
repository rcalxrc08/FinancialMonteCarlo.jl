# The following contains just abstract types regarding models and "virtual" method

#The most basic type, needed in order to support DifferentialEquations.jl
abstract type BaseProcess{T <: Number} end

# There will inherit from this type just processes, i.e. models that know what is a zeroCurve and a divided
abstract type AbstractMonteCarloProcess{T <: Number} <: BaseProcess{T} end

# Just engines will inherit from this, Brownian Motion et al.
abstract type AbstractMonteCarloEngine{T <: Number} <: BaseProcess{T} end

# Abstract type representing scalar montecarlo processes
abstract type ScalarMonteCarloProcess{T <: Number} <: AbstractMonteCarloProcess{T} end

# Abstract type representing multi dimensional montecarlo processes
abstract type VectorialMonteCarloProcess{T <: Number} <: AbstractMonteCarloProcess{T} end

# Abstract type representing N dimensional montecarlo processes (??)
abstract type NDimensionalMonteCarloProcess{T <: Number} <: VectorialMonteCarloProcess{T} end

# Abstract type for Levy (useful in FinancialFFT.jl)
abstract type LevyProcess{T <: Number} <: ScalarMonteCarloProcess{T} end

# Abstract type for Ito, an Ito process is always a Levy
abstract type ItoProcess{T <: Number} <: LevyProcess{T} end

# Abstract type for FA processes
abstract type FiniteActivityProcess{T <: Number} <: LevyProcess{T} end

# Abstract type for IA processes
abstract type InfiniteActivityProcess{T <: Number} <: LevyProcess{T} end

# Utility for dividends (implemented just for mono dimensional processes)
dividend(x::mc) where {mc <: ScalarMonteCarloProcess} = x.underlying.d;

"""
    simulate(mcProcess, zeroCurve, mcConfig, T)

Simulate the paths of the underlying of `mcProcess` from time 0 to `T`, under the risk neutral
measure given by `zeroCurve`.

Returns an `Nsim × (Nstep + 1)` matrix whose `i`-th row is the `i`-th path on the uniform time grid
`0:T/Nstep:T` (the first column is the spot). For multivariate processes, a vector with one such
matrix per component. The element type is inferred from the inputs, so it can be a dual number.

New models implement the in-place `simulate!(S, mcProcess, zeroCurve, mcConfig, T)`.
"""
function simulate(mcProcess::AbstractMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    price_type = predict_output_type_zero(mcProcess, zeroCurve, mcBaseData, T)
    S = get_matrix_type(mcBaseData, mcProcess, price_type)
    simulate!(S, mcProcess, zeroCurve, mcBaseData, T)
    return S
end

function simulate(mcProcess::VectorialMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    price_type = predict_output_type_zero(mcProcess, zeroCurve, mcBaseData, T)
    matrix_type = get_matrix_type(mcBaseData, mcProcess, price_type)
    S = matrix_type(undef, length(mcProcess.models))
    for i = 1:length(mcProcess.models)
        S[i] = eltype(S)(undef, mcBaseData.Nsim, mcBaseData.Nstep + 1)
    end
    simulate!(S, mcProcess, zeroCurve, mcBaseData, T)
    return S
end

"""
    simulate(engine, mcConfig, T)

Simulate the paths of an engine (e.g. [`BrownianMotion`](@ref)), which needs no zero rate curve,
from time 0 to `T`. Returns an `Nsim × (Nstep + 1)` matrix, one path per row.
"""
function simulate(mcProcess::AbstractMonteCarloEngine, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    price_type = predict_output_type_zero(mcProcess, mcBaseData, T)
    S = get_matrix_type(mcBaseData, mcProcess, price_type)
    simulate!(S, mcProcess, mcBaseData, T)
    return S
end

# function simulate_per_path(Nsim,Nstep,rn)
# S0=100.0;
# X=Matrix{Float64}(undef,Nsim,Nstep+1);
# for i=1:Nsim
# @views simulate_path!(X[i,:],rn)
# end

# return X;
# end

# function simulate_per_path(mcProcess::AbstractMonteCarloProcess, zeroCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
#     price_type = predict_output_type_zero(mcProcess, zeroCurve, mcBaseData, T)
#     S = get_matrix_type(mcBaseData, mcProcess, price_type)
#     #simulate!(S, mcProcess, zeroCurve, mcBaseData, T)
# 	for i=1:mcBaseData.Nsim
# 		# @views simulate_path!(X[i,:],rn)
# 		@views simulate_path!(S[i,:], mcProcess, zeroCurve, mcBaseData, T)
# 	end
#     return S
# end