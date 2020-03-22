using Statistics, DatagenCopulaBased 
"""
Struct for MultiVariate (log moneyness x-> S0*exp(x) ) Copula Process

		gaussianCopulaNVariateLogProcess=GaussianCopulaNVariateLogProcess(models::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		models  =	the processes.
		rho  = 	correlation matrix.
"""
mutable struct GaussianCopulaNVariateLogProcess{ num3 <: Number} <: NDimensionalMonteCarloProcess
	models::Tuple{Vararg{BaseProcess}}
	rho::Matrix{num3}
	function GaussianCopulaNVariateLogProcess(rho::Matrix{num3},models::BaseProcess...) where { num3 <: Number} 
		sz=size(rho)
		@assert sz[1]==sz[2]
		@assert length(models)==sz[1]
		@assert isposdef(rho)
		return new{num3}(models,rho);
	end
	function GaussianCopulaNVariateLogProcess(models::BaseProcess...) 
		len_=length(models)
		return GaussianCopulaNVariateLogProcess(Matrix{Float64}(I, len_, len_),models...);
	end
	function GaussianCopulaNVariateLogProcess(model1::BaseProcess,model2::BaseProcess,rho::num3) where { num3 <: Number} 
		corr_matrix_=[1.0 rho; rho 1.0];
		@assert det(corr_matrix_)>=0
		return GaussianCopulaNVariateLogProcess(corr_matrix_,model1,model2);
	end
end
 
export GaussianCopulaNVariateLogProcess;

function simulate!(S_Total,mcProcess::GaussianCopulaNVariateLogProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	@assert T>0.0
	
	####Simulation
	## Simulate
	len_=length(mcProcess.models);
	#S_Total::Array{Matrix{Number}}=[log.(simulate(model_i,rfCurve,mcBaseData,T)./model_i.underlying.S0) for model_i in mcProcess.models];
	
	for i in 1:len_
		simulate!(S_Total[i],mcProcess.models[i],rfCurve,mcBaseData,T)
	end
	
	rho=mcProcess.rho
	if (rho[1,2:end]==zeros(len_-1))
		return;
	end
	for i in 1:len_
		S_Total[i].=log.(S_Total[i]./mcProcess.models[i].underlying.S0);
	end
	for j in 1:Nstep
		U_joint=gausscopulagen2(Nsim,rho);
	
		for i in 1:len_
			cdf_ = sort(S_Total[i][:,j+1])
			@views S_Total[i][:,j+1]=(mcProcess.models[i].underlying.S0).*exp.(Statistics.quantile(cdf_,U_joint[:,i]))
		end
	end
	
	## Conclude
	return;

end
