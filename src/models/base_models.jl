# The following contains just abstract types regarding models and "virtual" method

abstract type BaseProcess end

abstract type AbstractMonteCarloProcess <: BaseProcess end

abstract type AbstractMonteCarloEngine <: BaseProcess end

abstract type ScalarMonteCarloProcess <: AbstractMonteCarloProcess end

abstract type VectorialMonteCarloProcess <: AbstractMonteCarloProcess end

abstract type NDimensionalMonteCarloProcess <: VectorialMonteCarloProcess end

abstract type LevyProcess<:ScalarMonteCarloProcess end

abstract type ItoProcess<:LevyProcess end

abstract type FiniteActivityProcess<:LevyProcess end

abstract type InfiniteActivityProcess<:LevyProcess end

dividend(x::mc) where {mc <: ScalarMonteCarloProcess} = x.underlying.d;

# export AbstractMonteCarloProcess
# export ScalarMonteCarloProcess
# export VectorialMonteCarloProcess
# export NDimensionalMonteCarloProcess
# export ItoProcess
# export LevyProcess
# export FiniteActivityProcess
# export InfiniteActivityProcess


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
function simulate(mcProcess::AbstractMonteCarloProcess,zeroCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	price_type=predict_output_type_zero(mcProcess,zeroCurve,mcBaseData,T);
	matrix_type=get_matrix_type(mcBaseData,mcProcess,price_type);
	S=matrix_type(undef,mcBaseData.Nsim,mcBaseData.Nstep+1);
	simulate!(S,mcProcess,zeroCurve,mcBaseData,T)
	return S;
end

function simulate(mcProcess::VectorialMonteCarloProcess,zeroCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	price_type=predict_output_type_zero(mcProcess,zeroCurve,mcBaseData,T);
	matrix_type=get_matrix_type(mcBaseData,mcProcess,price_type);
	S=matrix_type(undef,length(mcProcess.models));
	for i=1:length(mcProcess.models)
		S[i]=eltype(S)(undef,mcBaseData.Nsim,mcBaseData.Nstep+1);
	end
	simulate!(S,mcProcess,zeroCurve,mcBaseData,T)
	return S;
end

"""
General Interface for Stochastic Engine simulation

		S=simulate(mcProcess,mcBaseData,T)
	
Where:\n
		mcProcess          = Process to be simulated.
		mcBaseData = Basic properties of MonteCarlo simulation
		T = Final time of the process

		S      = Matrix with path of underlying.

"""
function simulate(mcProcess::AbstractMonteCarloEngine,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	price_type=predict_output_type_zero(mcProcess,mcBaseData,T);
	matrix_type=get_matrix_type(mcBaseData,mcProcess,price_type);
	S=matrix_type{typeof(price_type)}(undef,mcBaseData.Nsim,mcBaseData.Nstep+1);
	simulate!(S,mcProcess,mcBaseData,T)
	return S;
end