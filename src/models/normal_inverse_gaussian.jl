"""
Struct for NIG Process

		nigProcess=NormalInverseGaussianProcess(σ::num1,θ::num2,κ::num3) where {num1,num2,num3<: Number}
	
Where:\n
		σ =	volatility of the process.
		θ = variance of volatility.
		κ =	skewness of volatility.
"""
mutable struct NormalInverseGaussianProcess{num <: Number, num1 <: Number,num2<:Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: InfiniteActivityProcess{numtype}
	σ::num
	θ::num1
	κ::num2
	underlying::abstrUnderlying
	function NormalInverseGaussianProcess(σ::num,θ::num1,κ::num2,underlying::abstrUnderlying) where {num <: Number, num1 <: Number,num2 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ<=0
			error("volatility must be positive");
		elseif κ<=0
			error("κappa must be positive");
		elseif 1.0-(σ^2+2.0*θ)*κ<0
			error("Parameters with unfeasible values")
		else
			zero_typed=zero(num)+zero(num1)+zero(num2);
            return new{num,num1,num2,abstrUnderlying,typeof(zero_typed)}(σ,θ,κ,underlying)
        end
    end
end

export NormalInverseGaussianProcess;

function simulate!(X,mcProcess::NormalInverseGaussianProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nstep=mcBaseData.Nstep;

	σ=mcProcess.σ;
	θ1=mcProcess.θ;
	κ1=mcProcess.κ;
	@assert T>0
	
	dt=T/Nstep;
	psi1=(1-sqrt(1-(σ^2+2*θ1)*κ1))/κ1;
	drift=r-d-psi1;
	
	#Define subordinator
	IGRandomVariable=InverseGaussian(dt,dt*dt/κ1);
	
	#Call SubordinatedBrownianMotion
	simulate!(X,SubordinatedBrownianMotion(σ,drift,IGRandomVariable),mcBaseData,T);

	f(x)=S0*exp(x);
	broadcast!(f,X,X)

	return nothing;
	
end
