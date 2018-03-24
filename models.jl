abstract type AbstractMonteCarloProcess end


function simulate(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64)
	error("This is a Virtual function")
end

include("models/black_scholes.jl")
include("models/kou.jl")
include("models/merton.jl")
include("models/variance_gamma.jl")
include("models/normal_inverse_gaussian.jl")