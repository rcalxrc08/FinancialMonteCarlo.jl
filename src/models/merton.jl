type MertonProcess<:FiniteActivityProcess end

export MertonProcess;

function simulate(mcProcess::MertonProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,monteCarloMode::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	
	if(length(mcBaseData.param)!=4)
		error("Merton Model needs 4 parameters")
	end
	sigma=mcBaseData.param["sigma"];
	lambda1=mcBaseData.param["lambda"];
	mu1=mcBaseData.param["muJ"];
	sigma1=mcBaseData.param["sigmaJ"];
	if T<=0.0
		error("Final time must be positive");
	elseif lambda1<=0.0
		error("jump intensity must be positive");
	elseif sigma1<=0.0
		error("Final time must be positive");
	end
	
	####Simulation
	t=linspace(0,T,Nstep+1);
	## Simulate
	# -psi(-i)
	drift_RN=r-d-sigma^2/2-lambda1*(exp(mu1+sigma1*sigma1/2.0)-1.0); 
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift_RN)
	brownianMcData=MonteCarloBaseData(dict1,Nsim,Nstep);
	X=simulate(BrownianMotion(),spotData,brownianMcData,T,monteCarloMode)
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
