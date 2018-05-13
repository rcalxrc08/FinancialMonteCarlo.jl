type AmericanOption<:PayoffMC end

struct AMOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
end

export AmericanOption,AMOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,AmericanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.
		AmericanOption = Type of the Option
		isCall = true for Call Options, false for price Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},amOptionData::AMOptionData,spotData::equitySpotData,Payoff::AmericanOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	S0=S[1,1];
	Nsim=length(S[1:end,1]);
	Nstep=length(S[1,1:end])
	K=amOptionData.K;
	r=spotData.r;
	T=amOptionData.T;
	dt=T/Nstep
	# initialize 
	exerciseTimes=Nstep.*ones(Nsim);
	V=max.(0.0,iscall.*(S[:,end].-K)); #payoff
	# Backward Procedure 
	for j in Nstep-1:-1:1
		inMoneyIndexes=find((S[:,j].-K).*iscall.>0.0);
		if !isempty(inMoneyIndexes)
			S_I=S[inMoneyIndexes,j];
			#-- Intrinsic Value
			IV=(S_I-K)*iscall;
			#-- Continuation Value 
			#- Linear Regression on Quadratic Form
			A=[ones(length(S_I),1) S_I S_I.^2];
			b=V[inMoneyIndexes].*exp.(-r*dt*(exerciseTimes[inMoneyIndexes]-j));
			#MAT=A'*A;			
			LuMat=lufact(A'*A);
			Btilde=A'*b;
			alpha=LuMat\Btilde;
			#alpha=A\b;
			#Continuation Value
			CV=A*alpha;
			#----------
			# Find premature exercise times
			Index=find(IV.>CV);
			exercisePositions=inMoneyIndexes[Index];
			# Update the outputs
			V[exercisePositions]=IV[Index];
			exerciseTimes[exercisePositions]=j;
		end
	end
	price=max.(iscall*(S0-K),V.*exp.(-r*dt.*exerciseTimes))
	
	return price;
end
