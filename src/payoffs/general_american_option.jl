
struct AMOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
	function AMOptionData(T::Float64,K::Float64)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new(T,K)
        end
    end
end

export AMOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.
		isCall = true for Call Options, false for price Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},spotData::equitySpotData,phi::Function,T::Number) where{num<:Number}
	S0=S[1,1];
	(Nsim,Nstep)=size(S)
	Nstep-=1;
	r=spotData.r;
	dt=T/Nstep
	# initialize 
	exerciseTimes=Nstep.*ones(Nsim);
	#V=max.(0.0,iscall.*(S[:,end].-K)); #payoff
	V=phi.(S[:,end]); #payoff
	# Backward Procedure 
	for j in Nstep-1:-1:1
		inMoneyIndexes=find(phi.(S[:,j]).>0.0);
		if !isempty(inMoneyIndexes)
			S_I=S[inMoneyIndexes,j];
			#-- Intrinsic Value
			IV=phi.(S_I);
			#-- Continuation Value 
			#- Linear Regression on Quadratic Form
			A=[ones(length(S_I)) S_I S_I.^2];
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
	price=max.(phi(S0),V.*exp.(-r*dt.*exerciseTimes))
	
	return price;
end
