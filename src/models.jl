abstract type AbstractMonteCarloProcess end

export AbstractMonteCarloProcess

using Distributions;

function simulate(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64)
	error("This is a Virtual function")
end

include("models/brownian_motion.jl")
include("models/geometric_brownian_motion.jl")
include("models/black_scholes.jl")
include("models/kou.jl")
include("models/merton.jl")
include("models/subordinated_brownian_motion.jl")
include("models/variance_gamma.jl")
include("models/normal_inverse_gaussian.jl")