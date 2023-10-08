
abstract type AbstractMonteCarloConfiguration end

struct MonteCarloConfiguration{num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode} <: AbstractMonteCarloConfiguration
    Nsim::num1
    Nstep::num2
    monteCarloMethod::abstractMonteCarloMethod
    parallelMode::baseMode
    init_rng::Bool
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, seed::Number) where {num1 <: Integer, num2 <: Integer}
        return MonteCarloConfiguration(Nsim, Nstep, StandardMC(), SerialMode(seed))
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod = StandardMC(), seed::Number = 0) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod}
        return MonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, SerialMode(Int64(seed), MersenneTwister()))
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, parallelMethod::baseMode) where {num1 <: Integer, num2 <: Integer, baseMode <: BaseMode}
        return MonteCarloConfiguration(Nsim, Nstep, StandardMC(), parallelMethod)
    end
    #Most General, no default argument
    function GeneralMonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, init_rng::Bool = true) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode}
        ChainRulesCore.@ignore_derivatives @assert Nsim > zero(num1) "Number of Simulations must be positive"
        ChainRulesCore.@ignore_derivatives @assert Nstep > zero(num2) "Number of Steps must be positive"
        return new{num1, num2, abstractMonteCarloMethod, baseMode}(Nsim, Nstep, monteCarloMethod, parallelMethod, init_rng)
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, init_rng::Bool) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode}
        ChainRulesCore.@ignore_derivatives @assert Nsim > zero(num1) "Number of Simulations must be positive"
        ChainRulesCore.@ignore_derivatives @assert Nstep > zero(num2) "Number of Steps must be positive"
        return GeneralMonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, parallelMethod, init_rng)
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode}
        return GeneralMonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, parallelMethod)
    end
    function MonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::AntitheticMC, parallelMethod::baseMode = SerialMode()) where {num1 <: Integer, num2 <: Integer, baseMode <: BaseMode}
        ChainRulesCore.@ignore_derivatives @assert iseven(Nsim) "Antithetic support only even number of simulations"
        return GeneralMonteCarloConfiguration(Nsim, Nstep, monteCarloMethod, parallelMethod)
    end
end
SerialMonteCarloConfig = MonteCarloConfiguration{num1, num2, num3, ser_mode} where {num1 <: Integer, num2 <: Integer, num3 <: FinancialMonteCarlo.AbstractMonteCarloMethod, ser_mode <: FinancialMonteCarlo.SerialMode{rng_}} where {rng_ <: Random.AbstractRNG}
SerialAntitheticMonteCarloConfig = MonteCarloConfiguration{num1, num2, num3, ser_mode} where {num1 <: Integer, num2 <: Integer, num3 <: FinancialMonteCarlo.AntitheticMC, ser_mode <: FinancialMonteCarlo.SerialMode{rng_}} where {rng_ <: Random.AbstractRNG}
SerialSobolMonteCarloConfig = MonteCarloConfiguration{num1, num2, num3, ser_mode} where {num1 <: Integer, num2 <: Integer, num3 <: FinancialMonteCarlo.SobolMode, ser_mode <: FinancialMonteCarlo.SerialMode{rng_}} where {rng_ <: Random.AbstractRNG}

export MonteCarloConfiguration;
