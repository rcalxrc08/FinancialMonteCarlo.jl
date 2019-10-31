
function payoff(S::AbstractMatrix{num},spotData::equitySpotData,phi::Function,T::num2) where{num <: Number,num2 <: Number}
	S0=first(S);
	(Nsim,Nstep)=size(S)
	Nstep-=1;
	r=spotData.r;
	dt=T/Nstep
	# initialize 
	exerciseTimes=Nstep.*ones(Nsim);
	V=phi.(S[:,end]); #payoff
	# Backward Procedure 
	for j in Nstep:-1:2
		inMoneyIndexes=findall(phi.(S[:,j]).>0.0);
		if !isempty(inMoneyIndexes)
			S_I=S[inMoneyIndexes,j];
			#-- Intrinsic Value
			IV=phi.(S_I);
			#-- Continuation Value 
			#- Linear Regression on Quadratic Form
			A=[ones(length(S_I)) S_I S_I.^2];
			b=V[inMoneyIndexes].*exp.(-r*dt*(exerciseTimes[inMoneyIndexes].-j));
			#MAT=A'*A;			
			LuMat=lu(A'*A);
			Btilde=A'*b;
			alpha=LuMat\Btilde;
			#alpha=A\b;
			#Continuation Value
			CV=A*alpha;
			#@show size(CV)
			#----------
			# Find premature exercise times
			Index=findall(IV.>CV);
			exercisePositions=inMoneyIndexes[Index];
			# Update the outputs
			V[exercisePositions]=IV[Index];
			exerciseTimes[exercisePositions].=j-1;
		end
	end
	price=V.*exp.(-r*dt.*exerciseTimes)
	
	return price;
end
