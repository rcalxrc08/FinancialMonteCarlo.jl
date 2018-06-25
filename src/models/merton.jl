
type MertonProcess{num,num1,num2,num3<:Number}<:FiniteActivityProcess
	sigma::num
	lambda::num1
	muJ::num2
	sigmaJ::num3
	function MertonProcess(sigma::num,lambda::num1,muJ::num2,sigmaJ::num3) where {num,num1,num2,num3 <: Number}
        if sigma<=0.0
			error("volatility must be positive");
		elseif lambda<=0.0
			error("jump intensity must be positive");
		elseif sigmaJ<=0.0
			error("positive lambda must be positive");
		else
            return new{num,num1,num2,num3}(sigma,lambda,muJ,sigmaJ)
        end
    end
end

export MertonProcess;

function simulate(mcProcess::MertonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;

	sigma=mcProcess.sigma;
	lambda1=mcProcess.lambda;
	mu1=mcProcess.muJ;
	sigma1=mcProcess.sigmaJ;
	if T<=0.0
		error("Final time must be positive");
	end
	
	####Simulation
	t=linspace(0,T,Nstep+1);
	## Simulate
	# -psi(-i)
	drift_RN=r-d-sigma^2/2-lambda1*(exp(mu1+sigma1*sigma1/2.0)-1.0); 
	brownianMcData=MonteCarloConfiguration(Nsim,Nstep);
	X=simulate(BrownianMotion(sigma,drift_RN),spotData,brownianMcData,T,monteCarloMode)
	D1=Poisson(lambda1*T);
	NJumps=Int.(quantile.(D1,rand(Nsim)));
	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			#i=findfirst(x->x>=tjump,t);
			#jump_size=mu1+sigma1*randn()
			#X[ii,i:end]+=jump_size; #add jump component
			
			for i in 1:Nstep
			   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
				  jump_size=mu1+sigma1*randn() #Compute jump size
				  X[ii,i+1:end]+=jump_size; #add jump component
				  break;
			   end
			end
			
		end
	end
	## Conclude
	S=S0*exp.(X);
	return S;
	
end
