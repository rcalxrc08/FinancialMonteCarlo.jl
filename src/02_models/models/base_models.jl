# The following contains just abstract types regarding models and "virtual" method

#The most basic type, needed in order to support StochasticDiffEq.jl (MonteCarloDiffEqModel)
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

"""
    simulate(mcProcess, underlying, zeroCurve, mcConfig, T)

Simulate the paths of `underlying` (an [`Underlying`](@ref), a tuple of them for a multivariate model, one per
component) under `mcProcess` from time 0 to `T`, under the risk neutral measure given by `zeroCurve`.

Returns an `Nsim × (Nstep + 1)` matrix whose `i`-th row is the `i`-th path on the uniform time grid
`0:T/Nstep:T` (the first column is the spot). For multivariate processes, a vector with one such
matrix per component. The element type is inferred from the inputs, so it can be a dual number.

New models implement the in-place `simulate!(S, mcProcess, underlying, zeroCurve, mcConfig, T)`.
"""
function simulate(mcProcess::AbstractMonteCarloProcess, und::AbstractUnderlying, zeroCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    S = get_matrix_type(mcBaseData, mcProcess, price_type(mcProcess, und, zeroCurve, T))
    simulate!(S, mcProcess, und, zeroCurve, mcBaseData, T)
    return S
end

function simulate(mcProcess::VectorialMonteCarloProcess, unds::UnderlyingTuple, zeroCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    check_underlyings(mcProcess, unds)
    matrix_type = get_matrix_type(mcBaseData, mcProcess, price_type(mcProcess, unds, zeroCurve, T))
    S = matrix_type(undef, length(mcProcess.models))
    for i = 1:length(mcProcess.models)
        S[i] = eltype(S)(undef, mcBaseData.Nsim, mcBaseData.Nstep + 1)
    end
    simulate!(S, mcProcess, unds, zeroCurve, mcBaseData, T)
    return S
end

# One underlying per component of a multivariate model
check_underlyings(mcProcess::VectorialMonteCarloProcess, unds) = check_dimension(length(unds) == length(mcProcess.models), "one underlying per component of the model: $(length(mcProcess.models)) models, $(length(unds)) underlyings")
ChainRulesCore.@non_differentiable check_underlyings(::Any, ::Any)

"""
    simulate(engine, mcConfig, T)

Simulate the paths of an engine (e.g. [`BrownianMotion`](@ref)), which needs no zero rate curve,
from time 0 to `T`. Returns an `Nsim × (Nstep + 1)` matrix, one path per row.
"""
function simulate(mcProcess::AbstractMonteCarloEngine, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    S = get_matrix_type(mcBaseData, mcProcess, price_type(mcProcess, T))
    simulate!(S, mcProcess, mcBaseData, T)
    return S
end
