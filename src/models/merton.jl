type MertonProcess<:AbstractMonteCarloProcess end

export MertonProcess;

using Distributions;
function simulate(mcProcess::MertonProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64)
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
	####Simulation
	t=linspace(0,T,Nstep+1);
	## Simulate
	# -psi(-i)
	#calcolo drift in RN, cambia a seconda dell' AF che scelgo
	drift_RN=r-d+(-sigma^2/2-lambda1*(exp(mu1+sigma1*sigma1/2.0)-1.0)); 
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift_RN)
	brownianMcData=MonteCarloBaseData(dict1,Nsim,Nstep);
	X=simulate(BrownianMotion(),spotData,brownianMcData,T)
	D1=Poisson(lambda1*T);
	NJumps=Int.(quantile.(D1,rand(Nsim))); #invariante
	#Z=randn(Nsim,Nstep);
	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			#i=findfirst(x->x>=tjump,t);
			#jump_size=mu1+sigma1*randn()
			#X[ii,i:end]+=jump_size; #aggiungo componente di salto.
			
			for i in 1:Nstep
			   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
				  #u=rand(); #simulate Uniform([0,1])
				  jump_size=mu1+sigma1*randn() #Compute jump size
				  X[ii,i+1:end]+=jump_size; #aggiungo componente di salto.
				  break;
			   end
			end
			
		end
	end
	## Conclude
	S=S0*exp.(X);
	return S;
	
end
