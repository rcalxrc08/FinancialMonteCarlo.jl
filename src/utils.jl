
struct equitySpotData{T1,T2,T3<:Number}
    S0::T1
    r::T2
    d::T3
    function equitySpotData(S0::T1,r::T2,d::T3) where {T1, T2, T3 <: Number}
        if S0 <= 0.0
            error("Spot price must be positive")
        else
            return new{T1,T2,T3}(S0,r,d)
        end
    end
end 


abstract type BaseMode end
struct SerialMode <: BaseMode end
abstract type ParallelMode <: BaseMode end

abstract type AbstractMethod end
abstract type AbstractMonteCarloMethod <: AbstractMethod end
struct StandardMC <: AbstractMonteCarloMethod end
struct AntitheticMC <: AbstractMonteCarloMethod end

export equitySpotData,MonteCarloMode;

struct MonteCarloConfiguration{num1 <: Integer , num2 <: Integer , abstractMonteCarloMethod <: AbstractMonteCarloMethod }
	Nsim::num1
	Nstep::num2
	monteCarloMethod::abstractMonteCarloMethod
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,monteCarloMethod::abstractMonteCarloMethod=StandardMC()) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		else
            return new{num1,num2,abstractMonteCarloMethod}(Nsim,Nstep,monteCarloMethod)
        end
    end
end

export MonteCarloConfiguration;