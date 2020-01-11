abstract type BaseProcess end
abstract type AbstractUnderlying end

mutable struct Underlying{num <: Number, num2 <: Number} <: AbstractUnderlying
	S0::num
	d::num2
	function Underlying(S0::num_,d::num_2=0.0) where {num_ <: Number, num_2 <: Number}
		if(S0<num_(0))
			error("Underlying starting value must be positive");
		else
			return new{num_,num_2}(S0,d)
		end
	end
end

mutable struct UnderlyingVec{num <: Number, num2 <: Number, num3 <: Number} <: AbstractUnderlying
	S0::num
	d::Curve{num2,num3}
	function UnderlyingVec(S0::num_,d::Curve{num2,num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
		if(S0<num_(0))
			error("Underlying starting value must be positive");
		else
			return new{num_,num2,num3}(S0,d)
		end
	end
end

function Underlying(S0::num_,d::Curve{num2,num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
	return UnderlyingVec(S0,d)
end


export Underlying;
export UnderlyingVec;
	
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

using Distributions;

include("models/utils.jl")

### Ito Processes
include("models/brownian_motion.jl")
include("models/brownian_motion_aug.jl")
include("models/brownian_motion_prescribed.jl")
include("models/geometric_brownian_motion.jl")
include("models/geometric_brownian_motion_aug.jl")
include("models/black_scholes.jl")
include("models/heston.jl")
include("models/heston_aug.jl")
include("models/log_normal_mixture.jl")
include("models/shifted_log_normal_mixture.jl")

### Finite Activity Levy Processes
include("models/fa_levy.jl")
include("models/kou.jl")
include("models/merton.jl")

### Infinite Activity Levy Processes
include("models/subordinated_brownian_motion.jl")
include("models/subordinated_brownian_motion_aug.jl")
include("models/variance_gamma.jl")
include("models/normal_inverse_gaussian.jl")


include("models/nvariate.jl")
include("models/nvariate_log.jl")


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
function simulate(mcProcess::BaseProcess,zeroCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	error("Function used just for documentation")
end

include("models/operations.jl")