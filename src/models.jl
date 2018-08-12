import DiffEqBase.AbstractMonteCarloProblem

const BaseProcess=DiffEqBase.AbstractMonteCarloProblem

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


### Support for DiffentialEquations.jl
include("models/diff_eq_monte_carlo.jl")

############### Display Function

import Base.Multimedia.display;

function display(p::AbstractMonteCarloProcess)
	fldnames=fieldnames(typeof(p));
	for name in fldnames
		println(name," = ",getfield(p,name))
	end
end