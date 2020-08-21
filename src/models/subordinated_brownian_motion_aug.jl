"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1)
	
Where:\n
		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate,Continuous}, numtype <: Number} <: AbstractMonteCarloProcess{numtype}
	sigma::num
	drift::Curve{num1,num4}
	subordinator_::Distr
	function SubordinatedBrownianMotionVec(sigma::num,drift::Curve{num1,num4},dist::Distr) where {num <: Number,num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate,Continuous}}
        if sigma<=0
			error("volatility must be positive");
		else
			zero_typed=zero(num)+zero(num1)+zero(num4);
            return new{num,num1,num4,Distr,typeof(zero_typed)}(sigma,drift,dist)
        end
    end
end

function SubordinatedBrownianMotion(σ::num,drift::Curve{num1,num4},subordinator_::Distr) where {num <: Number, num1 <: Number, num4 <: Number,Distr <: Distribution{Univariate,Continuous}}
	return SubordinatedBrownianMotionVec(σ,drift,subordinator_)
end

export SubordinatedBrownianMotionVec;

function simulate!(X,mcProcess::SubordinatedBrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	@assert T>0
	type_sub=typeof(quantile(mcProcess.subordinator_,0.5));
	dt_s=Array{type_sub}(undef,Nsim)
	t_s=zeros(type_sub,Nsim)
	dt=T/Nstep;
	zero_drift=drift(zero(type_sub),zero(type_sub)+dt)*0;
	isDualZero=sigma*zero(type_sub)*0*zero_drift;
	@views X[:,1].=isDualZero;
	Z=Array{Float64}(undef,Nsim)
	for i=1:Nstep
		randn!(mcBaseData.rng,Z)
		rand!(mcBaseData.rng,mcProcess.subordinator_,dt_s);
		tmp_drift=[drift(t_,dt_) for (t_,dt_) in zip(t_s,dt_s)];
		t_s+=dt_s;
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		@views X[:,i+1]=X[:,i]+tmp_drift+sigma*sqrt.(dt_s).*Z;
	end

	return;
end


function simulate!(X,mcProcess::SubordinatedBrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	type_sub=typeof(quantile(mcProcess.subordinator_,0.5));
	dt_s=Array{type_sub}(undef,Nsim)
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	@assert T>0
	t_s=zeros(type_sub,Nsim)
	dt=T/Nstep;
	zero_drift=drift(zero(type_sub),zero(type_sub)+dt);
	isDualZero=sigma*zero(type_sub)*0*zero_drift;
	#X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	@views X[:,1].=isDualZero;
	for i=1:Nstep
		NsimAnti=div(Nsim,2)
		Z=randn(mcBaseData.rng,NsimAnti);
		#this is bad
		Z=[Z;-Z];
		rand!(mcBaseData.rng,mcProcess.subordinator_,dt_s);
		tmp_drift=[drift(t_,dt_) for (t_,dt_) in zip(t_s,dt_s)];
		t_s+=dt_s;
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		@views X[:,i+1]=X[:,i].+tmp_drift.+sigma.*sqrt.(dt_s).*Z;
	end
	
	return;
end
