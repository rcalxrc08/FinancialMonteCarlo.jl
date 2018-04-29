type HestonProcess<:ItoProcess end

export HestonProcess;

function simulate(mcProcess::HestonProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,mode1::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=6)
		error("Heston Model needs 6 parameters")
	end
	sigma=mcBaseData.param["sigma"];
	sigma_zero=mcBaseData.param["sigma_zero"];
	lambda1=mcBaseData.param["lambda"];
	kappa=mcBaseData.param["kappa"];
	rho=mcBaseData.param["rho"];
	theta=mcBaseData.param["theta"];

	####Simulation
	## Simulate
	kappa_s=kappa+lambda1;
	theta_s=kappa*theta/(kappa+lambda1);

	dt=T/Nstep
	isDualZero=S0*T*r*sigma_zero*kappa*theta*lambda1*sigma*rho*0.0;
	X=Matrix{typeof(isDualZero)}(Nsim,Nstep+1);
	X[:,1]=isDualZero;
	if mode1!=antithetic
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
