abstract type AbstractMonteCarloProcess end

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

function simulate(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,mode1::MonteCarloMode=standard)
	error("This is a Virtual function")
end

include("models/utils.jl")
include("models/brownian_motion.jl")
include("models/geometric_brownian_motion.jl")
include("models/black_scholes.jl")
include("models/kou.jl")
include("models/merton.jl")
include("models/subordinated_brownian_motion.jl")
include("models/variance_gamma.jl")
include("models/normal_inverse_gaussian.jl")
include("models/heston.jl")

include("models/diff_eq_monteCarlo.jl")