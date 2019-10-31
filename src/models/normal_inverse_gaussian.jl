"""
Struct for NIG Process

		nigProcess=NormalInverseGaussianProcess(σ::num1,θ::num2,κ::num3) where {num1,num2,num3<: Number}
	
Where:\n
		σ =	volatility of the process.
		θ = variance of volatility.
		κ =	skewness of volatility.
"""
mutable struct NormalInverseGaussianProcess{num <: Number, num1 <: Number,num2<:Number}<:InfiniteActivityProcess
	σ::num
	θ::num1
	κ::num2
	underlying::Underlying
	function NormalInverseGaussianProcess(σ::num,θ::num1,κ::num2,underlying::Underlying) where {num <: Number, num1 <: Number,num2 <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif κ<=0.0
			error("κappa must be positive");
		elseif 1.0-(σ^2+2.0*θ)*κ<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2}(σ,θ,κ,underlying)
        end
    end
	function NormalInverseGaussianProcess(σ::num,θ::num1,κ::num2,S0::num3) where {num <: Number, num1 <: Number,num2 <: Number, num3 <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif κ<=0.0
			error("κappa must be positive");
		elseif 1.0-(σ^2+2.0*θ)*κ<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2}(σ,θ,κ,Underlying(S0))
        end
    end
end

export NormalInverseGaussianProcess;

function simulate(mcProcess::NormalInverseGaussianProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	r=spotData.r;
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
	dt_s=quantile.([IGRandomVariable],rand(Nsim,Nstep));
	
	X=simulate(SubordinatedBrownianMotion(σ,drift,Underlying(0.0)),spotData,mcBaseData,T,dt_s);

	S=S0.*exp.(X);
	
	return S;
	
end
