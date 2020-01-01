using Distributed

import Future.randjump

abstract type AbstractZeroRateCurve end

struct ZeroRateCurve{T2 <: Number} <: AbstractZeroRateCurve
    r::T2
    function ZeroRateCurve(r::T2) where {T2 <: Number}
       return new{T2}(r)
    end
end
########### Curve Section

import Base.insert!
insert!(dic::Dict,k::Number,val::Number) = dic[k]=val; return nothing;
using Dictionaries
#const FinMCDict=Dict
const FinMCDict=HashDictionary
using Interpolations
function intgral_2(x::num,T::Array{num1},r::Array{num2}) where {num <: Number, num1 <: Number, num2 <: Number}
	if(x==0.0)
		return 0.0;
	end
	@assert x<=T[end]
	idx_=findfirst(y->y>=x,T)-1;
	out=sum([(r[i]+r[i+1])*0.5*(T[i+1]-T[i]) for i in 1:(idx_-1)])
	itp=LinearInterpolation([T[idx_],T[idx_+1]],[r[idx_],r[idx_+1]]);
	out=out+(r[idx_]+itp(x))*0.5*(x-T[idx_]);
	return out;
end

struct ZeroRateCurve2{num1 <: Number,num2 <: Number} <: AbstractZeroRateCurve
    r::FinMCDict{num1,num2}
	function ZeroRateCurve2(r_::Array{num1},T::num2) where {num1 <: Number, num2 <: Number}
		#@assert T==sort(T)
		r=FinMCDict{num2,num1}();
		Nstep=length(r_)-1;
		dt=T/Nstep;
		for i in 1:length(r_)
			#r[(i-1)*dt]=r_[i];
			insert!(r,(i-1)*dt,r_[i]);
		end
       new{num2,num1}(r)
    end
	function ZeroRateCurve2(r_::FinMCDict{num1,num2}) where {num1 <: Number, num2 <: Number}
       new{num2,num1}(r_)
    end
	function (x::ZeroRateCurve2)(t::Number,dt::Number)
		T=collect(keys(x.r));
		r=collect(values(x.r));
       return intgral_2(t+dt,T,r)-intgral_2(t,T,r);
    end
end 
function integral(x::ZeroRateCurve2,t::Number)
	T=collect(keys(x.r));
	r=collect(values(x.r));
   return intgral_2(t,T,r);
end
function integral(r::FinMCDict{num1,num2},t::Number) where {num1 <: Number, num2 <: Number}
	T=collect(keys(r));
	r=collect(values(r));
   return intgral_2(t,T,r);
end
integral(x::num1,t::num2) where {num1 <: Number, num2 <: Number}=x*t;

### end curve section

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


mutable struct MonteCarloConfiguration{num1 <: Integer , num2 <: Integer , abstractMonteCarloMethod <: AbstractMonteCarloMethod , baseMode <: BaseMode}
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