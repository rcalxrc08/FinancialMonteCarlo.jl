"""
Struct for NIG Process

		nigProcess=NormalInverseGaussianProcess(σ::num1,θ::num2,κ::num3) where {num1,num2,num3<: Number}
	
Where:\n
		σ =	volatility of the process.
		θ = variance of volatility.
		κ =	skewness of volatility.
"""
mutable struct NormalInverseGaussianProcess{num <: Number, num1 <: Number,num2<:Number, abstrUnderlying <: AbstractUnderlying}<:InfiniteActivityProcess
	σ::num
	θ::num1
	κ::num2
	underlying::abstrUnderlying
	function NormalInverseGaussianProcess(σ::num,θ::num1,κ::num2,underlying::abstrUnderlying) where {num <: Number, num1 <: Number,num2 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ<=0.0
			error("volatility must be positive");
		elseif κ<=0.0
			error("κappa must be positive");
		elseif 1.0-(σ^2+2.0*θ)*κ<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2,abstrUnderlying}(σ,θ,κ,underlying)
        end
    end
end

export NormalInverseGaussianProcess;

function simulate(mcProcess::NormalInverseGaussianProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;

	σ=mcProcess.σ;
	θ1=mcProcess.θ;
	κ1=mcProcess.κ;
	if T<=0.0
		error("Final time must be positive");
	end
	
	dt=T/Nstep;
	#psi1(v::Number)::Number=(1-sqrt(1.0+ ((σ*v)^2-2.0*1im*θ1*v)*κ1))/κ1;
	psi1=(1-sqrt(1.0-(σ^2+2.0*θ1)*κ1))/κ1;
	drift=r-d-psi1;
	
	#Simulate subordinator
	IGRandomVariable=InverseGaussian(dt,dt*dt/κ1);
	#dt_s=quantile.(IGRandomVariable,rand(mcBaseData.rng,Nsim,Nstep));
	
	X=simulate(SubordinatedBrownianMotion(σ,drift,IGRandomVariable,Underlying(0.0)),mcBaseData,T);

	S=S0.*exp.(X);
	
	return S;
	
end
