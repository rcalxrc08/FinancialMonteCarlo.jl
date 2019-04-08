using DifferentialEquations

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


@enum MonteCarloMode standard=1 antithetic=2


export equitySpotData,MonteCarloMode;

struct MonteCarloConfiguration{num1,num2<:Integer}
	Nsim::num1
	Nstep::num2
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2) where {num1,num2<:Integer}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		else
            return new{num1,num2}(Nsim,Nstep)
        end
    end
end

export MonteCarloConfiguration;