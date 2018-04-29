type KouProcess<:FiniteActivityProcess end

export KouProcess;

function simulate(mcProcess::KouProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,mode1::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=5)
		error("Kou Model needs 5 parameters")
	end
	sigma=mcBaseData.param["sigma"];
	lambda1=mcBaseData.param["lambda"];
	p=mcBaseData.param["p"];
	lambdap=mcBaseData.param["lambdap"];
	lambdam=mcBaseData.param["lambdam"];

	####Simulation
	## Simulate
	# r-d-psi(-i)
	drift_RN=r-d-sigma^2/2-lambda1*(p/(lambdap-1)-(1-p)/(lambdam+1));
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift_RN)
	brownianMcData=MonteCarloBaseData(dict1,Nsim,Nstep);
	X=simulate(BrownianMotion(),spotData,brownianMcData,T,mode1)

	t=linspace(0.0,T,Nstep+1);

	PoissonRV=Poisson(lambda1*T);
	PosExpRV=Exponential(1.0/lambdap);
	NegExpRV=Exponential(1.0/lambdam);
	NJumps=quantile.(PoissonRV,rand(Nsim)); #invariante

	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			
			idx1=findfirst(x->x>=tjump,t);
			u=rand(); #simulate Uniform([0,1])
			jump_size=(u<p)?quantile(PosExpRV,rand()):-quantile(NegExpRV,rand())
			X[ii,idx1:end]+=jump_size; #aggiungo componente di salto.
			
			#for i in 1:Nstep
			#   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
			#	  u=rand(); #simulate Uniform([0,1])
			#	  jump_size=(u<p)?quantile(PosExpRV,rand()):-quantile(NegExpRV,rand()) #Compute jump size
			#	  X[ii,i+1:end]+=jump_size; #aggiungo componente di salto.
			#	  break;
			#   end
			#end
			
		end
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end
