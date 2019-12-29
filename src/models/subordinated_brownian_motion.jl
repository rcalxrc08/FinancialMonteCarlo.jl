"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1,underlying::Underlying{nums0,numd})
	
Where:\n
		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotion{num <: Number, num1 <: Number, Distr <: Distribution{Univariate,Continuous}, nums0 <: Number, numd <: Number}<:AbstractMonteCarloProcess
	sigma::num
	drift::num1
	subordinator_::Distr
	underlying::Underlying{nums0,numd}
	function SubordinatedBrownianMotion(sigma::num,drift::num1,dist::Distr,underlying::Underlying{nums0,numd}) where {num <: Number,num1 <: Number, Distr <: Distribution{Univariate,Continuous}, nums0 <: Number, numd <: Number}
        if sigma<=0.0
			error("volatility must be positive");
		else
            return new{num,num1,Distr,nums0,numd}(sigma,drift,dist,underlying)
        end
    end
end

export SubordinatedBrownianMotion;

function simulate(mcProcess::SubordinatedBrownianMotion,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	dt_s=quantile.(mcProcess.subordinator_,rand(mcBaseData.rng,Nsim,Nstep));
	isDualZero=drift*sigma*dt_s[1,1]*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;

	Z=Array{Float64}(undef,Nsim)
	for i=1:Nstep
		randn!(mcBaseData.rng,Z)
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
	end

	return X;
end


function simulate(mcProcess::SubordinatedBrownianMotion,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	dt_s=quantile.(mcProcess.subordinator_,rand(mcBaseData.rng,Nsim,Nstep));
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	isDualZero=drift*sigma*dt_s[1,1]*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	for i=1:Nstep
		NsimAnti=div(Nsim,2)
		Z=randn(mcBaseData.rng,NsimAnti);
		Z=[Z;-Z];
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
	end
	
	return X;
end
