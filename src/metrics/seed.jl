import Future.randjump

randjump_(rng::MersenneTwister, num) = randjump(rng, num);
randjump_(rng, num) = rng;

inner_seed!(tmp_rng, seed) = Random.seed!(tmp_rng, seed)

function set_seed(mcConfig::MonteCarloConfiguration{type1, type2, type3, type4, type5}) where {type1 <: Integer, type2 <: Integer, type3 <: AbstractMonteCarloMethod, type4 <: SerialMode, type5 <: Random.AbstractRNG}
    tmp_rng = deepcopy(mcConfig.rng)
    inner_seed!(tmp_rng, mcConfig.seed)
    mcConfig.rng = randjump_(tmp_rng, mcConfig.offset)

    return
end
