
function payoff(S1::abstractMatrix,amPayoff::AmericanPayoff,spotData::ZeroRateCurve,T1::num2=amPayoff.T) where{abstractMatrix <: AbstractMatrix{num}} where{num <: Number, num2 <: Number}
	T=amPayoff.T;	
	(Nsim,NStep)=size(S1)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S=collect(S1[:,1:index1])
	S0=first(S);
	(Nsim,Nstep)=size(S)
	Nstep-=1;
	r=spotData.r;
	dt=T/Nstep
	# initialize 
	exerciseTimes=Nstep.*ones(Nsim);
	phi(Sti::numtype_) where {numtype_<:Number}=payout(Sti,amPayoff);
	@views V=phi.(S[:,end]); #payoff
	# Backward Procedure 
	for j in Nstep:-1:2
		@views Tmp=phi.(S[:,j]);
		inMoneyIndexes=findall(Tmp.>0.0);
		if !isempty(inMoneyIndexes)
			#S_I=S[inMoneyIndexes,j];
			S_I=view(S,inMoneyIndexes,j);
			#-- Intrinsic Value
			@views IV=Tmp[inMoneyIndexes];
			#-- Continuation Value 
			#- Linear Regression on Quadratic Form
			A=[ones(length(S_I)) S_I S_I.^2];
			@views b=V[inMoneyIndexes].*exp.(-r*dt*(exerciseTimes[inMoneyIndexes].-j));
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
			@views exercisePositions=inMoneyIndexes[Index];
			# Update the outputs
			@views V[exercisePositions]=IV[Index];
			@views exerciseTimes[exercisePositions].=j-1;
		end
	end
	Out=V.*exp.(-r*dt.*exerciseTimes)
	
	return Out;
end
