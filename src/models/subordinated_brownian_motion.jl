"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1,underlying::abstrUnderlying)
	
Where:\n
		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotion{num <: Number, num1 <: Number, Distr <: Distribution{Univariate,Continuous}, abstrUnderlying <: AbstractUnderlying}<:AbstractMonteCarloProcess
	sigma::num
	drift::num1
	subordinator_::Distr
	underlying::abstrUnderlying
	function SubordinatedBrownianMotion(sigma::num,drift::num1,dist::Distr,underlying::abstrUnderlying) where {num <: Number,num1 <: Number, Distr <: Distribution{Univariate,Continuous}, abstrUnderlying <: AbstractUnderlying}
        if sigma<=0.0
			error("volatility must be positive");
		else
            return new{num,num1,Distr,abstrUnderlying}(sigma,drift,dist,underlying)
        end
    end
end

export SubordinatedBrownianMotion;

function simulate(mcProcess::SubordinatedBrownianMotion,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	@assert T>0.0
	type_sub=typeof(quantile(mcProcess.subordinator_,0.5));
	dt_s=Array{type_sub}(undef,Nsim)
	isDualZero=drift*zero(type_sub)*0.0*mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	@views X[:,1].=isDualZero;
	Z=Array{Float64}(undef,Nsim)
	for i=1:Nstep
		randn!(mcBaseData.rng,Z)
		dt_s.=quantile.(mcProcess.subordinator_,rand(mcBaseData.rng,Nsim));
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		@views X[:,i+1].=X[:,i].+drift.*dt_s.+sigma.*sqrt.(dt_s).*Z;
	end

	return X;
end


function simulate(mcProcess::SubordinatedBrownianMotion,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	type_sub=typeof(quantile(mcProcess.subordinator_,0.5));
	dt_s=Array{type_sub}(undef,Nsim)
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	@assert T>0.0
	isDualZero=drift*sigma*zero(type_sub)*0.0*mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	for i=1:Nstep
		NsimAnti=div(Nsim,2)
		Z=randn(mcBaseData.rng,NsimAnti);
		Z=[Z;-Z];
		dt_s.=quantile.(mcProcess.subordinator_,rand(mcBaseData.rng,Nsim));
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i].+drift.*dt_s.+sigma.*sqrt.(dt_s).*Z;
	end
	
	return X;
end
