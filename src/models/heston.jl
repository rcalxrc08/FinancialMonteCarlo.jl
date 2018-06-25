
type HestonProcess{num,num1,num2,num3,num4,num5<:Number}<:ItoProcess
	sigma::num
	sigma_zero::num1
	lambda::num2
	kappa::num3
	rho::num4
	theta::num5
	function HestonProcess(sigma::num,sigma_zero::num1,lambda::num2,	kappa::num3,rho::num4,theta::num5) where {num,num1,num2,num3,num4,num5 <: Number}
        if sigma_zero<=0.0
			error("initial volatility must be positive");
		elseif sigma<=0.0
			error("volatility of volatility must be positive");
		elseif abs(kappa+lambda)<=1e-14
			error("unfeasible parameters");
		elseif !(-1.0<=rho<=1.0)
			error("rho must be a correlation");
        else
            return new{num,num1,num2,num3,num4,num5}(sigma,sigma_zero,lambda,kappa,rho,theta)
        end
    end
end

export HestonProcess;

function simulate(mcProcess::HestonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	sigma=mcProcess.sigma;
	sigma_zero=mcProcess.sigma_zero;
	lambda1=mcProcess.lambda;
	kappa=mcProcess.kappa;
	rho=mcProcess.rho;
	theta=mcProcess.theta;
	if T<=0.0
		error("Final time must be positive");
	end

	####Simulation
	## Simulate
	kappa_s=kappa+lambda1;
	theta_s=kappa*theta/(kappa+lambda1);

	dt=T/Nstep
	isDualZero=S0*T*r*sigma_zero*kappa*theta*lambda1*sigma*rho*0.0;
	X=Matrix{typeof(isDualZero)}(Nsim,Nstep+1);
	X[:,1]=isDualZero;
	if monteCarloMode!=antithetic
		v_m=[sigma_zero+isDualZero for _ in 1:Nsim];
		for j in 1:Nstep
			e1=randn(Nsim);
			e2=randn(Nsim);
			e2=e1.*rho.+e2.*sqrt(1-rho*rho);
			X[:,j+1]=X[:,j]+((r-d).-0.5.*max.(v_m,0)).*dt+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
			v_m+=kappa_s.*(theta_s.-max.(v_m,0)).*dt+sigma.*sqrt.(max.(v_m,0)).*sqrt(dt).*e2;
		end
	else
		v_m=[sigma_zero+isDualZero for _ in 1:Nsim];
		NsimAnti=Int(floor(Nsim/2))
		for j in 1:Nstep
			e1=randn(NsimAnti);
			e2=randn(NsimAnti);
			e1=[e1;-e1];
			e2=[e2;-e2];
			e2=e1.*rho.+e2.*sqrt(1-rho*rho);
			X[:,j+1]=X[:,j]+((r-d).-0.5.*max.(v_m,0)).*dt+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
			v_m+=kappa_s.*(theta_s.-max.(v_m,0)).*dt+sigma.*sqrt.(max.(v_m,0)).*sqrt(dt).*e2;
		end
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end
