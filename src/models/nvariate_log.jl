using EmpiricalCDFs, DatagenCopulaBased 
"""
Struct for Kou Process

		kouProcess=GaussianCopulaNVariateLogProcess(models::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		models  =	volatility of the process.
		λ  = 	jumps intensity.
		p  =	prob. of having positive jump.
		λ₊ =	positive jump size.
		λ₋ =	negative jump size.
"""
mutable struct GaussianCopulaNVariateLogProcess{ num3 <: Number} <: NDimensionalMonteCarloProcess
	models::Tuple{Vararg{FinancialMonteCarlo.BaseProcess}}
	rho::Matrix{num3}
	function GaussianCopulaNVariateLogProcess(rho::Matrix{num3},models::FinancialMonteCarlo.BaseProcess...) where { num3 <: Number} 
		sz=size(rho)
		@assert sz[1]==sz[2]
		@assert length(models)==sz[1]
		return new{num3}(models,rho);
	end
end

export GaussianCopulaNVariateLogProcess;

function simulate(mcProcess::GaussianCopulaNVariateLogProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if T<=0.0
		error("Final time must be positive");
	end
	
	####Simulation
	## Simulate
	S_1=log.(simulate(mcProcess.models[1],rfCurve,mcBaseData,T)./mcProcess.models[1].underlying.S0);
	
	#Preallocate in some way
	S_Total=[S_1];
	len_=length(mcProcess.models);
	for i in 2:len_
		S_tmp=log.(simulate(mcProcess.models[i],rfCurve,mcBaseData,T)./mcProcess.models[i].underlying.S0);
		push!(S_Total,S_tmp)
	end
	
	rho=mcProcess.rho
	for j in 1:Nstep
	
		U_joint=gausscopulagen(Nsim,rho);
	
		for i in 1:len_
			S_tmp=S_Total[i][:,j+1];
			cdf_ = EmpiricalCDF()
			append!(cdf_,S_tmp)
			sort!(cdf_)
			icdf_ = finv(cdf_)
			S_Total[i][:,j+1]=(mcProcess.models[i].underlying.S0).*exp.(icdf_.(U_joint[:,i]))
		end
	end
	
	## Conclude
	return S_Total;

end
