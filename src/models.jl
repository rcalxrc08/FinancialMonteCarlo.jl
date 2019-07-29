abstract type BaseProcess end

abstract type AbstractMonteCarloProcess <: BaseProcess end

abstract type ItoProcess<:AbstractMonteCarloProcess end

abstract type LevyProcess<:AbstractMonteCarloProcess end

abstract type FiniteActivityProcess<:LevyProcess end

abstract type InfiniteActivityProcess<:LevyProcess end

export AbstractMonteCarloProcess
export ItoProcess
export LevyProcess
export FiniteActivityProcess
export InfiniteActivityProcess

using Distributions;

include("models/utils.jl")

### Ito Processes
include("models/brownian_motion.jl")
include("models/geometric_brownian_motion.jl")
include("models/black_scholes.jl")
include("models/heston.jl")
include("models/log_normal_mixture.jl")

### Finite Activity Levy Processes
include("models/kou.jl")
include("models/merton.jl")

### Infinite Activity Levy Processes
include("models/subordinated_brownian_motion.jl")
include("models/variance_gamma.jl")
include("models/normal_inverse_gaussian.jl")

############### Display Function

import Base.Multimedia.display;

function display(p::Union{AbstractMonteCarloProcess,AbstractPayoff})
	fldnames=collect(fieldnames(typeof(p)));
	print(typeof(p),"(");
	print(fldnames[1],"=",getfield(p,fldnames[1]))
	popfirst!(fldnames)
	for name in fldnames
		print(",",name,"=",getfield(p,name))
	end
	println(")");
end

"""
General Interface for Stochastic Process simulation

		S=simulate(mcProcess,spotData,mcBaseData,T,monteCarloMode=standard,parallelMode=SerialMode())
	
Where:\n
		mcProcess          = Process to be simulated.
		spotData  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		parallelMode  [Optional, default to SerialMode()] = SerialMode(), CudaMode(), AFMode()

		S      = Matrix with path of underlying.

"""
function simulate(mcProcess::BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3},T::numb,parallelMode::BaseMode=SerialMode()) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod}
	error("Function used just for documentation")
end