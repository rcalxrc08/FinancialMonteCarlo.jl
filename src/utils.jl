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

@enum MonteCarloMode standard=1 antithetic=2

export equitySpotData,MonteCarloMode;

struct MonteCarloConfiguration
	Nsim::Integer
	Nstep::Integer
	function MonteCarloConfiguration(Nsim::Integer,Nstep::Integer)
        if Nsim <= 0
            error("Number of Simulations must be positive")
        elseif Nstep <= 0
            error("Number of Steps must be positive")
		else
            return new(Nsim,Nstep)
        end
    end
end

export MonteCarloConfiguration;