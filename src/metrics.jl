function set_seed(mcConfig::MonteCarloConfiguration{type1,type2,type3,type4})  where {type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: SerialMode}

	Random.seed!(mcConfig.rng,mcConfig.seed)
	#mcConfig.rng=randjump_(mcConfig.rng, mcConfig.offset)
	
	return;
	
end


include("metrics/pricer.jl")
include("metrics/pricer_cv.jl")
include("metrics/delta.jl")
include("metrics/confinterval.jl")
include("metrics/variance.jl")