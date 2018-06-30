
type NormalInverseGaussianProcess{num,num1,num2<:Number}<:InfiniteActivityProcess
	σ::num
	θ::num1
	κ::num2
	function NormalInverseGaussianProcess(σ::num,θ::num1,κ::num2) where {num,num1,num2 <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif κ<=0.0
			error("κappa must be positive");
		elseif 1.0-(σ^2+2.0*θ)*κ<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2}(σ,θ,κ)
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
	dt_s=quantile.(IGRandomVariable,rand(Nsim,Nstep));
	
	X=simulate(SubordinatedBrownianMotion(σ,drift),spotData,mcBaseData,T,dt_s,monteCarloMode);

	S=S0.*exp.(X);
	
	return S;
	
end
