"""
Struct for Kou Process

		kouProcess=KouProcess(σ::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		p  =	prob. of having positive jump.
		λ₊ =	positive jump size.
		λ₋ =	negative jump size.
"""
mutable struct KouProcess{num <: Number, num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number}<:FiniteActivityProcess
	σ::num
	λ::num1
	p::num2
	λ₊::num3
	λ₋::num4
	function KouProcess(σ::num,λ::num1,p::num2,λ₊::num3,λ₋::num4) where {num <: Number, num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif λ<=0.0
			error("jump intensity must be positive");
		elseif λ₊<=0.0
			error("positive λ must be positive");
		elseif λ₋<=0.0
			error("negative λ must be positive");
		elseif !(0<=p<=1)
			error("p must be a probability")
        else
            return new{num,num1,num2,num3,num4}(σ,λ,p,λ₊,λ₋)
        end
    end
end

export KouProcess;

function simulate(mcProcess::KouProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	λ1=mcProcess.λ;
	p=mcProcess.p;
	λ₊=mcProcess.λ₊;
	λ₋=mcProcess.λ₋;
	if T<=0.0
		error("Final time must be positive");
	end

	####Simulation
	## Simulate
	# r-d-psi(-i)
	drift_RN=r-d-σ^2/2-λ1*(p/(λ₊-1)-(1-p)/(λ₋+1));
	X=Matrix(simulate(BrownianMotion(σ,drift_RN),spotData,mcBaseData,T,parallelMode))

	t=range(0.0, stop=T, length=Nstep+1);
	PoissonRV=Poisson(λ1*T);
	NJumps=quantile.([PoissonRV],rand(Nsim));

	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the number of jumps (conditional simulation)
		tjumps=sort(rand(Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			
			idx1=findfirst(x->x>=tjump,t);
			u=rand(); #simulate Uniform([0,1])
			jump_size=u<p ? quantile_exp(λ₊,rand()) : -quantile_exp(λ₋,rand())
			X[ii,idx1:end].+=jump_size; #add jump component
			
			#for i in 1:Nstep
			#   if tjump>t[i] && tjump<=t[i+1] #Look for where it is happening the jump
			#	  u=rand(); #simulate Uniform([0,1])
			#	  jump_size=(u<p) ? quantile(PosExpRV,rand()):-quantile(NegExpRV,rand()) #Compute jump size
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
