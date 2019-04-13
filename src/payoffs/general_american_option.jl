
function payoff(S::AbstractMatrix{num},spotData::equitySpotData,phi::Function,T::num2) where{num,num2<:Number}
	S0=spotData.S0;
	(Nsim,Nstep)=size(S)
	Nstep-=1;
	r=spotData.r;
	d=spotData.d;
	dt=T/Nstep
	# initialize 
	exerciseTimes=Nstep.*ones(Nsim);
	#V=max.(0.0,iscall.*(S[:,end].-K)); #payoff
	V=phi.(S[:,end]); #payoff
	# Backward Procedure 
	for j in Nstep-1:-1:1
		inMoneyIndexes=findall(phi.(S[:,j]).>0.0);
		if !isempty(inMoneyIndexes)
			S_I=S[inMoneyIndexes,j];
			#-- Intrinsic Value
			IV=phi.(S_I);
			#-- Continuation Value 
			#- Linear Regression on Quadratic Form
			A=[ones(length(S_I)) S_I S_I.^2];
			b=V[inMoneyIndexes].*exp.(-(r-d)*dt*(exerciseTimes[inMoneyIndexes].-j));
			#MAT=A'*A;			
			LuMat=lu(A'*A);
			Btilde=A'*b;
			alpha=LuMat\Btilde;
			#alpha=A\b;
			#Continuation Value
			CV=A*alpha;
			#----------
			# Find premature exercise times
			Index=findall(IV.>CV);
			exercisePositions=inMoneyIndexes[Index];
			# Update the outputs
			V[exercisePositions]=IV[Index];
			exerciseTimes[exercisePositions].=j;
		end
	end
	price=max.(phi(S0),V.*exp.(-d*dt.*exerciseTimes))
	
	return price;
end
