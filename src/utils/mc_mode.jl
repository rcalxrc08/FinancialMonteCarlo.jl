#Parallel Modes
abstract type BaseMode end

"""
    SerialMode(seed = 0, rng = MersenneTwister())

Single threaded simulation drawing from `rng`, which is reset to `seed` by every metric.
Any `Random.AbstractRNG` can be used.
"""
struct SerialMode{rngType <: Random.AbstractRNG} <: BaseMode
    seed::Int64
    rng::rngType
    function SerialMode(seed::num = 0, rng::rngType = MersenneTwister()) where {num <: Integer, rngType <: Random.AbstractRNG}
        return new{rngType}(Int64(seed), rng)
    end
end
abstract type ParallelMode <: BaseMode end

#Numerical Modes, keyholes for modes different than Monte Carlo
abstract type AbstractMethod end
abstract type AbstractMonteCarloMethod <: AbstractMethod end

"Plain Monte Carlo sampling (the default method of [`MonteCarloConfiguration`](@ref))."
struct StandardMC <: AbstractMonteCarloMethod end

"Antithetic variates: the second half of the paths uses the negated normals of the first half. Requires an even number of paths."
struct AntitheticMC <: AbstractMonteCarloMethod end

"Quasi Monte Carlo sampling with Sobol low discrepancy sequences."
struct SobolMode <: AbstractMonteCarloMethod end

"Reserved for simulations from prescribed random numbers; not implemented by the built-in models."
struct PrescribedMC <: AbstractMonteCarloMethod end
