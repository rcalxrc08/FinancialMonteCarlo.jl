
type KouProcess{num,num1,num2,num3,num4<:Number}<:FiniteActivityProcess
	sigma::num
	lambda::num1
	p::num2
	lambdap::num3
	lambdam::num4
	function KouProcess(sigma::num,lambda::num1,p::num2,lambdap::num3,lambdam::num4) where {num,num1,num2,num3,num4 <: Number}
        if sigma<=0.0
			error("volatility must be positive");
		elseif lambda<=0.0
			error("jump intensity must be positive");
		elseif lambdap<=0.0
			error("positive lambda must be positive");
		elseif lambdam<=0.0
			error("negative lambda must be positive");
		elseif !(0<=p<=1)
			error("p must be a probability")
        else
            return new{num,num1,num2,num3,num4}(sigma,lambda,p,lambdap,lambdam)
        end
    end
end

export KouProcess;

function simulate(mcProcess::KouProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	sigma=mcProcess.sigma;
	lambda1=mcProcess.lambda;
	p=mcProcess.p;
	lambdap=mcProcess.lambdap;
	lambdam=mcProcess.lambdam;
	if T<=0.0
		error("Final time must be positive");
	end

	####Simulation
	## Simulate
	# r-d-psi(-i)
	drift_RN=r-d-sigma^2/2-lambda1*(p/(lambdap-1)-(1-p)/(lambdam+1));
	brownianMcData=MonteCarloConfiguration(Nsim,Nstep);
	X=simulate(BrownianMotion(sigma,drift_RN),spotData,brownianMcData,T,monteCarloMode)

	t=linspace(0.0,T,Nstep+1);

	PoissonRV=Poisson(lambda1*T);
	NJumps=quantile.(PoissonRV,rand(Nsim));

	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			
			idx1=findfirst(x->x>=tjump,t);
			u=rand(); #simulate Uniform([0,1])
			jump_size=u<p?quantile_exp(lambdap,rand()):-quantile_exp(lambdam,rand())
			X[ii,idx1:end]+=jump_size; #add jump component
			
			#for i in 1:Nstep
			#   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
			#	  u=rand(); #simulate Uniform([0,1])
			#	  jump_size=(u<p)?quantile(PosExpRV,rand()):-quantile(NegExpRV,rand()) #Compute jump size
			#	  X[ii,i+1:end]+=jump_size; #add jump component
			#	  break;
			#   end
			#end
			
		end
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end
