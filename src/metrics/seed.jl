inner_seed!(rng, seed) = Random.seed!(rng, seed)

function set_seed!(mcConfig::MonteCarloConfiguration{type1, type2, type3, type4}) where {type1 <: Integer, type2 <: Integer, type3 <: AbstractMonteCarloMethod, type4 <: SerialMode}
    if mcConfig.init_rng
        inner_seed!(mcConfig.parallelMode.rng, mcConfig.parallelMode.seed)
    end
end
