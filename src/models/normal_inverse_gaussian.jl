
type NormalInverseGaussianProcess{num,num1,num2<:Number}<:InfiniteActivityProcess
	sigma::num
	theta::num1
	k::num2
	function NormalInverseGaussianProcess(sigma::num,theta::num1,k::num2) where {num,num1,num2 <: Number}
        if sigma<=0.0
			error("volatility must be positive");
		elseif k<=0.0
			error("kappa must be positive");
		elseif 1.0-(sigma^2+2.0*theta)*k<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2}(sigma,theta,k)
        end
    end
end

export NormalInverseGaussianProcess;

function simulate(mcProcess::NormalInverseGaussianProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;

	sigma=mcProcess.sigma;
	theta1=mcProcess.theta;
	k1=mcProcess.k;
	if T<=0.0
		error("Final time must be positive");
	end
	
	dt=T/Nstep;
	#psi1(v::Number)::Number=(1-sqrt(1.0+ ((sigma*v)^2-2.0*1im*theta1*v)*k1))/k1;
	psi1=(1-sqrt(1.0-(sigma^2+2.0*theta1)*k1))/k1;
	drift=r-d-psi1;
	subParam=MonteCarloConfiguration(Nsim,Nstep);
	
	#Simulate subordinator
	IGRandomVariable=InverseGaussian(dt,dt*dt/k1);
	dt_s=quantile.(IGRandomVariable,rand(Nsim,Nstep));
	
	X=simulate(SubordinatedBrownianMotion(sigma,drift),spotData,subParam,T,dt_s,monteCarloMode);

	S=S0.*exp.(X);
	
	return S;
	
end
