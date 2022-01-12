import Future.randjump

inner_seed!(tmp_rng, seed) = Random.seed!(tmp_rng, seed)

function set_seed(mcConfig::MonteCarloConfiguration{type1, type2, type3, type4}) where {type1 <: Integer, type2 <: Integer, type3 <: AbstractMonteCarloMethod, type4 <: SerialMode}
    tmp_rng = deepcopy(mcConfig.parallelMode.rng)
    inner_seed!(tmp_rng, mcConfig.parallelMode.seed)
    mcConfig.parallelMode.rng = tmp_rng

    return
end
