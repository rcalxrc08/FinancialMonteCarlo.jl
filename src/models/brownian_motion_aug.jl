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
	if σ <= 0.0
		error("Volatility must be positive")
	else
		return BrownianMotionVec(σ,μ,underlying)
	end
end

function simulate(mcProcess::BrownianMotionVec,rfCurve::ZeroRateCurve2,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	μ=mcProcess.μ;
	if T<=0.0
		error("Final time must be positive");
	end
	dt=T/Nstep
	stddev_bm=σ*sqrt(dt)
	isDualZero=stddev_bm*0.0+mcProcess.underlying.S0;
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