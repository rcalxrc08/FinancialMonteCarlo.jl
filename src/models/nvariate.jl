using EmpiricalCDFs, DatagenCopulaBased,LinearAlgebra
"""
Struct for MultiVariate Copula Process

		kouProcess=GaussianCopulaNVariateProcess(models::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		models  =	the processes.
		rho  = 	correlation matrix.
"""
mutable struct GaussianCopulaNVariateProcess{ num3 <: Number} <: NDimensionalMonteCarloProcess
	models::Tuple{Vararg{FinancialMonteCarlo.BaseProcess}}
	rho::Matrix{num3}
	function GaussianCopulaNVariateProcess(rho::Matrix{num3},models::FinancialMonteCarlo.BaseProcess...) where { num3 <: Number} 
		sz=size(rho)
		@assert sz[1]==sz[2]
		@assert length(models)==sz[1]
		@assert det(rho)>=0
		return new{num3}(models,rho);
	end
	function GaussianCopulaNVariateProcess(model1::FinancialMonteCarlo.BaseProcess,model2::FinancialMonteCarlo.BaseProcess,rho::num3) where { num3 <: Number} 
		corr_matrix_=[1.0 rho; rho 1.0];
		@assert det(corr_matrix_)>=0
		return new{num3}((model1,model2),[1.0 rho; rho 1.0]);
	end
end

export GaussianCopulaNVariateProcess;

function simulate(mcProcess::GaussianCopulaNVariateProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if T<=0.0
		error("Final time must be positive");
	end
	
	####Simulation
	## Simulate
	S_1=simulate(mcProcess.models[1],rfCurve,mcBaseData,T);
	
	#Preallocate in some way
	S_Total=[S_1];
	len_=length(mcProcess.models);
	for i in 2:len_
		S_tmp=simulate(mcProcess.models[i],rfCurve,mcBaseData,T);
		push!(S_Total,S_tmp)
	end
	
	rho=mcProcess.rho
	if (rho[1,2:end]==zeros(len_-1))
		return S_Total
	end
	for j in 1:Nstep
	
		U_joint=gausscopulagen(Nsim,rho);
	
		for i in 1:len_
			S_tmp=S_Total[i][:,j+1];
			cdf_ = EmpiricalCDF()
			append!(cdf_,S_tmp)
			sort!(cdf_)
			icdf_ = finv(cdf_)
			S_Total[i][:,j+1]=icdf_.(U_joint[:,i])
		end
	end
	
	## Conclude
	return S_Total;

end
