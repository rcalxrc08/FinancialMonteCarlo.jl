
abstract type AbstractMonteCarloConfiguration end

mutable struct MonteCarloConfiguration{num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode, rngType <: Random.AbstractRNG} <: AbstractMonteCarloConfiguration
    Nsim::num1
    Nstep::num2
    monteCarloMethod::abstractMonteCarloMethod
    parallelMode::baseMode
    seed::Int64
    offset::Union{Int64, Nothing}
    rng::rngType
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, seed::Number) where {num1 <: Integer, num2 <: Integer}
        return MonteCarloConfiguration(Nsim, Nstep, StandardMC(), SerialMode(), Int64(seed), MersenneTwister())
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod = StandardMC(), seed::Number = 0) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod}
        return MonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, SerialMode(), Int64(seed), MersenneTwister())
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, parallelMethod::baseMode, seed::Number = 0, rng::rngType_ = MersenneTwister()) where {num1 <: Integer, num2 <: Integer, baseMode <: BaseMode, rngType_ <: Random.AbstractRNG}
        return MonteCarloConfiguration(Nsim, Nstep, StandardMC(), parallelMethod, Int64(seed), rng)
    end
    #Most General, no default argument, offset is controllable from outside
    function GeneralMonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, seed::Number, rng::rngType_) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode, rngType_ <: Random.AbstractRNG}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
        else
            return new{num1, num2, abstractMonteCarloMethod, baseMode, rngType_}(Nsim, Nstep, monteCarloMethod, parallelMethod, Int64(seed), 0, rng)
        end
    end
    function GeneralMonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, seed::Number, rng::MersenneTwister, offset::Number) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
        else
            return new{num1, num2, abstractMonteCarloMethod, baseMode, MersenneTwister}(Nsim, Nstep, monteCarloMethod, parallelMethod, Int64(seed), offset, rng)
        end
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, seed::Number, rng::MersenneTwister, offset::Number) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
        else
            return GeneralMonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, parallelMethod, Int64(seed), offset, rng)
        end
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, seed::Number = 0, rng::rngType_ = MersenneTwister()) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode, rngType_ <: Random.AbstractRNG}
        return GeneralMonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, parallelMethod, Int64(seed), rng)
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::AntitheticMC, parallelMethod::baseMode, seed::Number = 0, rng::rngType_ = MersenneTwister()) where {num1 <: Integer, num2 <: Integer, baseMode <: BaseMode, rngType_ <: Random.AbstractRNG}
        if div(Nsim, 2) * 2 != Nsim
            error("Antithetic support only even number of simulations")
        else
            return GeneralMonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, parallelMethod, Int64(seed), rng)
        end
    end
end

export MonteCarloConfiguration;
