#Parallel Modes
abstract type BaseMode end
mutable struct SerialMode{rngType <: Random.AbstractRNG} <: BaseMode
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
struct StandardMC <: AbstractMonteCarloMethod end
struct AntitheticMC <: AbstractMonteCarloMethod end
struct SobolMode <: AbstractMonteCarloMethod end
struct PrescribedMC <: AbstractMonteCarloMethod end
