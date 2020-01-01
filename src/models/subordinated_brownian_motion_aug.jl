"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1,underlying::abstrUnderlying)
	
Where:\n
		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate,Continuous},  abstrUnderlying <: AbstractUnderlying}<:AbstractMonteCarloProcess
	sigma::num
	drift::Curve{num1,num4}
	subordinator_::Distr
	underlying::abstrUnderlying
	function SubordinatedBrownianMotionVec(sigma::num,drift::Curve{num1,num4},dist::Distr,underlying::abstrUnderlying) where {num <: Number,num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate,Continuous}, abstrUnderlying <: AbstractUnderlying}
        if sigma<=0.0
			error("volatility must be positive");
		else
            return new{num,num1,num4,Distr,abstrUnderlying}(sigma,drift,dist,underlying)
        end
    end
end

function SubordinatedBrownianMotion(σ::num,drift::Curve{num1,num4},subordinator_::Distr,underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num4 <: Number,Distr <: Distribution{Univariate,Continuous}, abstrUnderlying <: AbstractUnderlying}
	if σ <= 0.0
		error("Volatility must be positive")
	else
		return SubordinatedBrownianMotionVec(σ,drift,subordinator_,underlying)
	end
end

export SubordinatedBrownianMotionVec;

function simulate(mcProcess::SubordinatedBrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	type_sub=typeof(quantile(mcProcess.subordinator_,rand()));
	dt_s=Array{type_sub}(undef,Nsim)
	dt=T/Nstep;
	isDualZero=sigma*dt_s[1,1]*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	@views X[:,1].=isDualZero;
	Z=Array{Float64}(undef,Nsim)
	for i=1:Nstep
		randn!(mcBaseData.rng,Z)
		dt_s.=quantile.(mcProcess.subordinator_,rand(mcBaseData.rng,Nsim));
		tmp_drift=[drift((i-1)*dt,dt_) for dt_ in dt_s];
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		@views X[:,i+1].=X[:,i].+tmp_drift.*dt_s.+sigma.*sqrt.(dt_s).*Z;
	end

	return X;
end


function simulate(mcProcess::SubordinatedBrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	type_sub=typeof(quantile(mcProcess.subordinator_,rand()));
	dt_s=Array{type_sub}(undef,Nsim)
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	isDualZero=sigma*dt_s[1,1]*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	for i=1:Nstep
		NsimAnti=div(Nsim,2)
		Z=randn(mcBaseData.rng,NsimAnti);
		Z=[Z;-Z];
		dt_s.=quantile.(mcProcess.subordinator_,rand(mcBaseData.rng,Nsim));
		#I am not sure about this
		tmp_drift=[drift((i-1)*dt,dt_) for dt_ in dt_s];
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i].+tmp_drift.*dt_s.+sigma.*sqrt.(dt_s).*Z;
	end
	
	return X;
end
