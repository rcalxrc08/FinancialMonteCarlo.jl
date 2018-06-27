"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amPayoff)
	
Where:\n
		S           = Paths of the Underlying.
		amPayoff  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{num},spotData::equitySpotData,phi::Function,T::Number) where{num<:Number}
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
