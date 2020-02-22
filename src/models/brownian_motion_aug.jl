"""
Struct for Brownian Motion

		bmProcess=BrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct BrownianMotionVec{num <: Number, num1 <: Number , num4 <: Number, abstrUnderlying <: AbstractUnderlying} <: ItoProcess
	σ::num
	μ::Curve{num1,num4}
	underlying::abstrUnderlying
	function BrownianMotionVec(σ::num,μ::Curve{num1,num4},underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,num4,abstrUnderlying}(σ,μ,underlying)
        end
    end
end

function BrownianMotion(σ::num,μ::Curve{num1,num4},underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
	return BrownianMotionVec(σ,μ,underlying)
end

function simulate(mcProcess::BrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	μ=mcProcess.μ;
	@assert T>0.0
	dt=T/Nstep
	stddev_bm=σ*sqrt(dt)
	zero_drift=μ(dt*0.0,dt);
	isDualZero=stddev_bm*0.0*mcProcess.underlying.S0*zero_drift;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	view(X,:,1).=isDualZero;
	@inbounds for j=1:Nstep
		tmp_=μ((j-1)*dt,dt);
		@inbounds for i=1:Nsim
			x_i_j=@views X[i,j];
			@views X[i,j+1]=x_i_j+tmp_+stddev_bm*randn(mcBaseData.rng);
		end
	end

	return X;

end


function simulate(mcProcess::BrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	μ=mcProcess.μ;
	@assert T>0.0
	dt=T/Nstep
	stddev_bm=σ*sqrt(dt)
	zero_drift=μ(dt*0.0,dt);
	isDualZero=stddev_bm*0.0*mcProcess.underlying.S0*zero_drift;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	view(X,:,1).=isDualZero;
	Nsim_2=div(Nsim,2)

	@inbounds for j in 1:Nstep
		tmp_=μ((j-1)*dt,dt);
		@inbounds for i in 1:Nsim_2
			Z=stddev_bm*randn(mcBaseData.rng);
			X[2*i-1,j+1]=X[2*i-1,j]+tmp_+Z;
			X[2*i,j+1]  =X[2*i,j]  +tmp_-Z;
		end
	end

	return X;

end