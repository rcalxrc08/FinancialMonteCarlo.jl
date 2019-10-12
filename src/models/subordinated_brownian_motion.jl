
mutable struct SubordinatedBrownianMotion{num <: Number, num1 <: Number}<:AbstractMonteCarloProcess
	sigma::num
	drift::num1
	function SubordinatedBrownianMotion(sigma::num,drift::num1) where {num <: Number,num1 <: Number}
        if sigma<=0.0
			error("volatility must be positive");
		else
            return new{num,num1}(sigma,drift)
        end
    end
end

export SubordinatedBrownianMotion;

function simulate(mcProcess::SubordinatedBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3},T::numb,dt_s::Array{num1,2},parallelMode::SerialMode=SerialMode()) where {numb <: Number, num1<:Number , type1 <: Number, type2<: Number, type3 <: StandardMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(size(dt_s)!=(Nsim,Nstep))
		error("Inconsistent Time Matrix")
	end
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	isDualZero=drift*sigma*dt_s[1,1]*0.0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;

	Z=Array{Float64}(undef,Nsim)
	for i=1:Nstep
		randn!(Z)
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
	end

	return X;
end


function simulate(mcProcess::SubordinatedBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3},T::numb,dt_s::Array{num1,2},parallelMode::SerialMode=SerialMode()) where {numb <: Number, num1<:Number , type1 <: Number, type2<: Number, type3 <: AntitheticMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(size(dt_s)!=(Nsim,Nstep))
		error("Inconsistent Time Matrix")
	end
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	isDualZero=drift*sigma*dt_s[1,1]*0.0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	for i=1:Nstep
		NsimAnti=Int(floor(Nsim/2))
		Z=randn(NsimAnti);
		Z=[Z;-Z];
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
	end
	
	return X;
end
