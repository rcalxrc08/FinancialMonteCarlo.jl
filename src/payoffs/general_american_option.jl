mutable struct GeneralAmericanOption{amPayoff <: AmericanPayoff}<:AmericanPayoff
	option::amPayoff
	function GeneralAmericanOption(option::amPayoff) where {amPayoff <: AmericanPayoff}
		return new{amPayoff}(option)
    end
end


function payoff(S::abstractMatrix,amPayoff::GeneralAmericanOption,spotData::ZeroRateCurve,T::num2) where{abstractMatrix <: AbstractMatrix{num}} where{num <: Number, num2 <: Number}
	S0=first(S);
	(Nsim,Nstep)=size(S)
	Nstep-=1;
	r=spotData.r;
	dt=T/Nstep
	# initialize 
	exerciseTimes=Nstep.*ones(Nsim);
	phi(Sti::numtype_) where {numtype_<:Number}=payout(Sti,amPayoff.option);
	V=phi.(S[:,end]); #payoff
	# Backward Procedure 
	for j in Nstep:-1:2
		Tmp=phi.(S[:,j]);
		inMoneyIndexes=findall(Tmp.>0.0);
		if !isempty(inMoneyIndexes)
			S_I=S[inMoneyIndexes,j];
			#-- Intrinsic Value
			IV=Tmp[inMoneyIndexes];
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
			#----------
			# Find premature exercise times
			Index=findall(IV.>CV);
			exercisePositions=inMoneyIndexes[Index];
			# Update the outputs
			V[exercisePositions]=IV[Index];
			exerciseTimes[exercisePositions].=j-1;
		end
	end
	Out=V.*exp.(-r*dt.*exerciseTimes)
	
	return Out;
end
