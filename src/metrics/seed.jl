inner_seed!(rng, seed) = Random.seed!(rng, seed)

reseed!(state::SeededRNG) = inner_seed!(state.rng, state.seed)
reseed!(states::AbstractVector{<:SeededRNG}) = foreach(reseed!, states)

# The generators that a mode owns besides those of the configuration: the CUDA extension seeds the one of the device.
seed_device!(::AbstractMode, ::Integer) = nothing

function set_seed!(mcConfig::MonteCarloConfiguration)
    if mcConfig.reseed
        reseed!(mcConfig.rng)
        seed_device!(mcConfig.mode, base_seed(mcConfig.rng))
    end
    return nothing
end
