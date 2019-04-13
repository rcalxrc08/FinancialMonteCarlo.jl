if(DIFFEQ_MONTECARLO_ACTIVE_FLAG)
	import DiffEqBase.AbstractMonteCarloProblem
	const BaseProcess=DiffEqBase.AbstractMonteCarloProblem
else
	abstract type BaseProcess end
end

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

using Distributions,DualNumbers;

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


if(DIFFEQ_MONTECARLO_ACTIVE_FLAG)
	### Support for DiffentialEquations.jl
	include("models/diff_eq_monte_carlo.jl")
end

############### Display Function

import Base.Multimedia.display;

function display(p::AbstractMonteCarloProcess)
	fldnames=fieldnames(typeof(p));
	for name in fldnames
		println(name," = ",getfield(p,name))
	end
end

"""
General Interface for Stochastic Process simulation

		S=simulate(mcProcess,spotData,mcBaseData,T,monteCarloMode=standard,parallelMode=SerialMode())
	
Where:\n
		mcProcess          = Process to be simulated.
		spotData  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		monteCarloMode [Optional, default to standard]= standard or antitethic
		parallelMode  [Optional, default to SerialMode()] = SerialMode(), CudaMode(), AFMode()

		S      = Matrix with path of underlying.
```
"""
function simulate(mcProcess::BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode()) where {numb<:Number}
	error("Function used just for documentation")
end