# The following contains just abstract types regarding models and "virtual" method

abstract type BaseProcess end

abstract type AbstractMonteCarloProcess <: BaseProcess end

abstract type ScalarMonteCarloProcess <: AbstractMonteCarloProcess end

abstract type VectorialMonteCarloProcess <: AbstractMonteCarloProcess end

abstract type NDimensionalMonteCarloProcess <: VectorialMonteCarloProcess end

abstract type LevyProcess<:ScalarMonteCarloProcess end

abstract type ItoProcess<:LevyProcess end

abstract type FiniteActivityProcess<:LevyProcess end

abstract type InfiniteActivityProcess<:LevyProcess end

dividend(x::mc) where {mc <: ScalarMonteCarloProcess} = x.underlying.d;

export AbstractMonteCarloProcess
export ScalarMonteCarloProcess
export VectorialMonteCarloProcess
export NDimensionalMonteCarloProcess
export ItoProcess
export LevyProcess
export FiniteActivityProcess
export InfiniteActivityProcess


"""
General Interface for Stochastic Process simulation

		S=simulate(mcProcess,zeroCurve,mcBaseData,T)
	
Where:\n
		mcProcess          = Process to be simulated.
		zeroCurve  = Datas of the Zero Rate Curve.
		mcBaseData = Basic properties of MonteCarlo simulation
		T = Final time of the process

		S      = Matrix with path of underlying.

"""
function simulate(mcProcess::BaseProcess,zeroCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	error("Function used just for documentation")
end