function set_seed(mcConfig::MonteCarloConfiguration{type1,type2,type3,type4})  where {type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: SerialMode}

	tmp_rng=deepcopy(mcConfig.rng);
	Random.seed!(tmp_rng,mcConfig.seed)
	mcConfig.rng=randjump_(tmp_rng, mcConfig.offset)
	Random.seed!(mcConfig.seed*10+10)
	
	return;
	
end


include("metrics/pricer.jl")
include("metrics/pricer_cv.jl")
include("metrics/distribution.jl")
include("metrics/delta.jl")
include("metrics/rho.jl")
include("metrics/confinterval.jl")
include("metrics/variance.jl")