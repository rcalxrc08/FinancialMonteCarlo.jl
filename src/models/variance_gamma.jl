
type VarianceGammaProcess{num,num1,num2<:Number}<:InfiniteActivityProcess
	sigma::num
	theta::num1
	k::num2
	function VarianceGammaProcess(sigma::num,theta::num1,k::num2) where {num,num1,num2 <: Number}
        if sigma<=0.0
			error("volatility must be positive");
		elseif k<=0.0
			error("kappa must be positive");
		elseif 1-sigma*sigma*k/2.0-theta*k<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2}(sigma,theta,k)
        end
    end
end


export VarianceGammaProcess;

function simulate(mcProcess::VarianceGammaProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
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
	#-1/p[3]*log(1+u*u*p[1]*p[1]*p[3]/2.0-1im*p[2]*p[3]*u);
	psi1=-1/k1*log(1-sigma*sigma*k1/2.0-theta1*k1);
	#1-sigma*sigma*k1/2.0-theta1*k1
	drift=r-d-psi1;
	
	gammaRandomVariable=Gamma(dt/k1);
	dt_s=k1.*quantile.(gammaRandomVariable,rand(Nsim,Nstep));
	subParam=MonteCarloConfiguration(Nsim,Nstep);
	
	X=simulate(SubordinatedBrownianMotion(sigma,drift),spotData,subParam,T,dt_s,monteCarloMode);

	S=S0.*exp.(X);
	
	return S;
	
end
