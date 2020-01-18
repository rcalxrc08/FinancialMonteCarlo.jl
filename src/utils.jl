using Distributed

import Future.randjump

abstract type AbstractZeroRateCurve end

struct ZeroRate{T2 <: Number} <: AbstractZeroRateCurve
    r::T2
    function ZeroRate(r::T2) where {T2 <: Number}
       return new{T2}(r)
    end
end
########### Curve Section
import Base.insert!
insert!(dic::Dict,k::Number,val::Number) = dic[k]=val; return nothing;
using Dictionaries
#const FinMCDict=Dict
const FinMCDict=HashDictionary

const Curve = FinMCDict{num1,num2} where {num1 <: Number,num2 <: Number}

function keys_(x::Curve{num1,num2}) where {num1 <: Number,num2 <: Number}
	return sort(collect(keys(x)));
end

function values_(x::Curve{num1,num2}) where {num1 <: Number,num2 <: Number}
	T_d=collect(keys(x));
	idx_=sortperm(T_d);
	return collect(values(x))[idx_];
end

function Curve(r_::Array{num1},T::num2) where {num1 <: Number, num2 <: Number}
	r=FinMCDict{num2,num1}();
	Nstep=length(r_)-1;
	dt=T/Nstep;
	for i in 1:length(r_)
		insert!(r,(i-1)*dt,r_[i]);
	end
   return r
end
function Curve(r_::Array{num1},T::Array{num2}) where {num1 <: Number, num2 <: Number}
	@assert length(r_)==length(T)
	@assert T==sort(T)
	r=FinMCDict{num2,num1}();
	Nstep=length(r_)-1;
	for i in 1:length(r_)
		insert!(r,T[i],r_[i]);
	end
   return r
end
function ImpliedCurve(r_::Array{num1},T::Array{num2}) where {num1 <: Number, num2 <: Number}
	@assert length(r_)==length(T)
	@assert length(r_)>=1
	@assert T==sort(T)
	
	r=FinMCDict{num2,num1}();
	Nstep=length(r_)-1;
	insert!(r,0.0,0.0);
	insert!(r,T[1],r_[1]*2.0);
	prec_r=r_[1]*2.0;
	for i in 2:length(r_)
		tmp_r=(r_[i]*T[i]-r_[i-1]*T[i-1])*2/(T[i]-T[i-1])-prec_r
		insert!(r,T[i],tmp_r);
		prec_r=tmp_r;
	end
	return r;
end
function ImpliedCurve(r_::Array{num1},T::num2) where {num1 <: Number, num2 <: Number}
	N=length(r_);
	dt_=T/(N);
	t_=collect(dt_:dt_:T)
	return ImpliedCurve(r_,t_)
end
function (x::Curve)(t::Number,dt::Number)
	T=collect(keys_(x));
	r=collect(values_(x));
	return intgral_2(t+dt,T,r)-intgral_2(t,T,r);
end



using Interpolations
function intgral_2(x::num,T::Array{num1},r::Array{num2}) where {num <: Number, num1 <: Number, num2 <: Number}
	@assert length(T)==length(r)
	if(x==0.0)
		return 0.0;
	end
	idx_=findlast(y->y<x,T);
	out=sum([(r[i]+r[i+1])*0.5*(T[i+1]-T[i]) for i in 1:(idx_-1)])
	if x<=T[end]
		itp=LinearInterpolation([T[idx_],T[idx_+1]],[r[idx_],r[idx_+1]], extrapolation_bc = Flat());
		out=out+(r[idx_]+itp(x))*0.5*(x-T[idx_]);
	else
		#continuation
		out=out+(r[idx_]+r[idx_-1])*0.5*(x-T[idx_]);
	end
	return out;
end

struct ZeroRateCurve{num1 <: Number,num2 <: Number} <: AbstractZeroRateCurve
    r::Curve{num1,num2}
	function ZeroRateCurve(r_::Array{num1},T::num2) where {num1 <: Number, num2 <: Number}
       new{num2,num1}(Curve(r_,T))
    end

	function ZeroRateCurve(r_::Curve{num1,num2}) where {num1 <: Number, num2 <: Number}
       new{num2,num1}(r_)
    end
end 
function ZeroRate(r_::Array{num1},T::num2) where {num1 <: Number, num2 <: Number}
       return ZeroRateCurve(r_,T);
end
function ImpliedZeroRate(r_::Array{num1},T::num2) where {num1 <: Number, num2 <: Number}
   return ZeroRateCurve(ImpliedCurve(r_,T))
end
function ImpliedZeroRate(r_::Array{num1},T::Array{num2}) where {num1 <: Number, num2 <: Number}
   return ZeroRateCurve(ImpliedCurve(r_,T))
end

export ZeroRateCurve;

function integral(r::FinMCDict{num1,num2},t::Number) where {num1 <: Number, num2 <: Number}
	T=collect(keys_(r));
	r=collect(values_(r));
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
struct SobolMode <: AbstractMonteCarloMethod end
struct PrescribedMC <: AbstractMonteCarloMethod end

export ZeroRate;

function randjump_(rng,num)
	return randjump(rng,num);
end	


mutable struct MonteCarloConfiguration{num1 <: Integer , num2 <: Integer , abstractMonteCarloMethod <: AbstractMonteCarloMethod , baseMode <: BaseMode , rngType <: Random.AbstractRNG}
	Nsim::num1
	Nstep::num2
	monteCarloMethod::abstractMonteCarloMethod
	parallelMode::baseMode
	seed::Int64
	offset::Int64
	rng::rngType
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,seed::Number) where {num1 <: Integer, num2 <: Integer}
        return MonteCarloConfiguration(Nsim,Nstep,StandardMC(),SerialMode(),Int64(seed),MersenneTwister())
    end
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,monteCarloMethod::abstractMonteCarloMethod=StandardMC(),seed::Number=0) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod}
        return MonteCarloConfiguration(Nsim,Nstep,monteCarloMethod,SerialMode(),Int64(seed),MersenneTwister())
    end
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,parallelMethod::baseMode,seed::Number=0,rng::rngType_=MersenneTwister()) where {num1 <: Integer, num2 <: Integer, baseMode <: BaseMode, rngType_ <: Random.AbstractRNG}
        return MonteCarloConfiguration(Nsim,Nstep,StandardMC(),parallelMethod,Int64(seed),rng)
    end
	#Most General
	function MonteCarloConfiguration(Nsim::num1,Nstep::num2,monteCarloMethod::abstractMonteCarloMethod,parallelMethod::baseMode,seed::Number=0,rng::rngType_=MersenneTwister()) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: BaseMode, rngType_ <: Random.AbstractRNG}
        if Nsim <= zero(num1)
            error("Number of Simulations must be positive")
        elseif Nstep <= zero(num2)
            error("Number of Steps must be positive")
		elseif (abstractMonteCarloMethod <: AntitheticMC) & ( div(Nsim,2)*2!=Nsim )
			error("Antithetic support only even number of simulations")
		else
            return new{num1,num2,abstractMonteCarloMethod,baseMode,rngType_}(Nsim,Nstep,monteCarloMethod,parallelMethod,Int64(seed),div(Nsim,2)*Nstep*(Distributed.myid() == 1 ? 0 : (Distributed.myid()-2) ),rng)
        end
    end
end

export MonteCarloConfiguration;