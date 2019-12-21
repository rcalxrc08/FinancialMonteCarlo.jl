using Distributed

import Future.randjump
struct ZeroRateCurve{T2 <: Number}
    r::T2
    function ZeroRateCurve(r::T2) where {T2 <: Number}
       return new{T2}(r)
    end
end 


abstract type BaseMode end
struct SerialMode <: BaseMode end
abstract type ParallelMode <: BaseMode end

abstract type AbstractMethod end
abstract type AbstractMonteCarloMethod <: AbstractMethod end
struct StandardMC <: AbstractMonteCarloMethod end
struct AntitheticMC <: AbstractMonteCarloMethod end
struct PrescribedMC <: AbstractMonteCarloMethod end

export ZeroRateCurve;

function randjump_(rng,num)
	return randjump(rng,num);
end	


struct MonteCarloConfiguration{num1 <: Integer , num2 <: Integer , abstractMonteCarloMethod <: AbstractMonteCarloMethod , baseMode <: BaseMode}
	Nsim::num1
	Nstep::num2
	monteCarloMethod::abstractMonteCarloMethod
	parallelMode::baseMode
	seed::Int64
	offset::Int64
	rng::MersenneTwister
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,seed::Int64) where {num1 <: Integer, num2 <: Integer}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		else
            return new{num1,num2,StandardMC,SerialMode}(Nsim,Nstep,StandardMC(),SerialMode(),seed,div(Nsim,2)*Nstep*(Distributed.myid() == 1 ? 0 : (Distributed.myid()-2) ),MersenneTwister(seed))
        end
    end
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,monteCarloMethod::abstractMonteCarloMethod=StandardMC(),seed::Int64=0) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		elseif (abstractMonteCarloMethod <: AntitheticMC) & ( div(Nsim,2)*2!=Nsim )
			error("Antithetic support only even number of simulations")
		else
            return new{num1,num2,abstractMonteCarloMethod,SerialMode}(Nsim,Nstep,monteCarloMethod,SerialMode(),seed,div(Nsim,2)*Nstep*(Distributed.myid() == 1 ? 0 : (Distributed.myid()-2) ),MersenneTwister(seed))
        end
    end
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,parallelMethod::baseMode,seed::Int64=0) where {num1 <: Integer, num2 <: Integer, baseMode <: BaseMode}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		else
            return new{num1,num2,StandardMC,baseMode}(Nsim,Nstep,StandardMC(),parallelMethod,seed,div(Nsim,2)*Nstep*(Distributed.myid() == 1 ? 0 : (Distributed.myid()-2) ),MersenneTwister(seed))
        end
    end
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,monteCarloMethod::abstractMonteCarloMethod,parallelMethod::baseMode,seed::Int64=0) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		elseif (abstractMonteCarloMethod <: AntitheticMC) & ( div(Nsim,2)*2!=Nsim )
			error("Antithetic support only even number of simulations")
		else
            return new{num1,num2,abstractMonteCarloMethod,baseMode}(Nsim,Nstep,monteCarloMethod,parallelMethod,seed,div(Nsim,2)*Nstep*(Distributed.myid() == 1 ? 0 : (Distributed.myid()-2) ),MersenneTwister(seed))
        end
    end
end

export MonteCarloConfiguration;