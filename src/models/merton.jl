"""
Struct for Merton Process

		mertonProcess=MertonProcess(σ::num1,λ::num2,μJ::num3,σJ::num4) where {num1,num2,num3,num4<: Number}
	
Where:\n
		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		μJ =	jumps mean.
		σJ =	jumps variance.
"""
mutable struct MertonProcess{num <: Number, num1 <: Number,num2 <: Number,num3<:Number, nums0 <: Number, numd <: Number}<:FiniteActivityProcess
	σ::num
	λ::num1
	μJ::num2
	σJ::num3
	underlying::Underlying{nums0,numd}
	function MertonProcess(σ::num,λ::num1,μJ::num2,σJ::num3,underlying::Underlying{nums0,numd}) where {num <: Number, num1 <: Number,num2 <: Number,num3 <: Number, nums0 <: Number, numd <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif λ<=0.0
			error("jump intensity must be positive");
		elseif σJ<=0.0
			error("positive λ must be positive");
		else
            return new{num,num1,num2,num3,nums0,numd}(σ,λ,μJ,σJ,underlying)
        end
    end
	function MertonProcess(σ::num,λ::num1,μJ::num2,σJ::num3,S0::num4) where {num <: Number, num1 <: Number,num2 <: Number,num3 <: Number, num4 <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif λ<=0.0
			error("jump intensity must be positive");
		elseif σJ<=0.0
			error("positive λ must be positive");
		else
            return new{num,num1,num2,num3,num4,Float64}(σ,λ,μJ,σJ,Underlying(S0))
        end
    end
end

export MertonProcess;

function simulate(mcProcess::MertonProcess,spotData::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	r=spotData.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
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
	X=Matrix(simulate(BrownianMotion(σ,drift_RN,Underlying(0.0)),spotData,mcBaseData,T))
	D1=Poisson(λ1*T);
	NJumps=Int.(quantile.([D1],rand(mcBaseData.rng,Nsim)));
	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(mcBaseData.rng,Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			#i=findfirst(x->x>=tjump,t);
			#jump_size=μ+σ1*randn()
			#X[ii,i:end]+=jump_size; #add jump component
			
			for i in 1:Nstep
			   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
				  jump_size=μ+σ1*randn(mcBaseData.rng) #Compute jump size
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
