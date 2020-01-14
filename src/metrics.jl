function set_seed(mcConfig::MonteCarloConfiguration{type1,type2,type3,type4,type5})  where {type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: SerialMode, type5 <: Random.AbstractRNG}

	tmp_rng=deepcopy(mcConfig.rng);
	Random.seed!(tmp_rng,mcConfig.seed)
	if type5 <: MersenneTwister
		mcConfig.rng=randjump_(tmp_rng, mcConfig.offset)
	else
		mcConfig.rng=tmp_rng
	end
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