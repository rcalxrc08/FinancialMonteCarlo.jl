
mutable struct MertonProcess{num,num1,num2,num3<:Number}<:FiniteActivityProcess
	σ::num
	λ::num1
	μJ::num2
	σJ::num3
	function MertonProcess(σ::num,λ::num1,μJ::num2,σJ::num3) where {num,num1,num2,num3 <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif λ<=0.0
			error("jump intensity must be positive");
		elseif σJ<=0.0
			error("positive λ must be positive");
		else
            return new{num,num1,num2,num3}(σ,λ,μJ,σJ)
        end
    end
end

export MertonProcess;

function simulate(mcProcess::MertonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard) where {numb<:Number}
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;

	σ=mcProcess.σ;
	λ1=mcProcess.λ;
	μ=mcProcess.μJ;
	σ1=mcProcess.σJ;
	if T<=0.0
		error("Final time must be positive");
	end
	
	####Simulation
	t=range(0.0, stop=T, length=Nstep+1);
	## Simulate
	# -psi(-i)
	drift_RN=r-d-σ^2/2-λ1*(exp(μ+σ1*σ1/2.0)-1.0); 
	X=simulate(BrownianMotion(σ,drift_RN),spotData,mcBaseData,T,monteCarloMode)
	D1=Poisson(λ1*T);
	NJumps=Int.(quantile.(D1,rand(Nsim)));
	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			#i=findfirst(x->x>=tjump,t);
			#jump_size=μ+σ1*randn()
			#X[ii,i:end]+=jump_size; #add jump component
			
			for i in 1:Nstep
			   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
				  jump_size=μ+σ1*randn() #Compute jump size
				  X[ii,i+1:end].+=jump_size; #add jump component
				  break;
			   end
			end
			
		end
	end
	## Conclude
	S=S0*exp.(X);
	return S;
	
end
