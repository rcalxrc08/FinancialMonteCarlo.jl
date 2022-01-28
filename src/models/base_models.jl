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
General Interface for Stochastic Process simulation

		S=simulate(mcProcess,zeroCurve,mcBaseData,T)
	
Where:

		mcProcess   = Process to be simulated.
		zeroCurve   = Datas of the Zero Rate Curve.
		mcBaseData  = Basic properties of MonteCarlo simulation
		T           = Final time of the process

		S           = Matrix with path of underlying.

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
General Interface for Stochastic Engine simulation

		S=simulate(mcProcess,mcBaseData,T)
	
Where:

		mcProcess  = Process to be simulated.
		mcBaseData = Basic properties of MonteCarlo simulation
		T          = Final time of the process

		S          = Matrix with path of underlying.

"""
function simulate(mcProcess::AbstractMonteCarloEngine, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    price_type = predict_output_type_zero(mcProcess, mcBaseData, T)
    S = get_matrix_type(mcBaseData, mcProcess, price_type)
    simulate!(S, mcProcess, mcBaseData, T)
    return S
end
