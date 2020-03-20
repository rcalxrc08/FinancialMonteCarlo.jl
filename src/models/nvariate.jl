using Statistics, DatagenCopulaBased,LinearAlgebra
"""
Struct for MultiVariate Copula Process

		gaussianCopulaNVariateProcess=GaussianCopulaNVariateProcess(models::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		models  =	the processes.
		rho  = 	correlation matrix.
"""
mutable struct GaussianCopulaNVariateProcess{ num3 <: Number} <: NDimensionalMonteCarloProcess
	models::Tuple{Vararg{BaseProcess}}
	rho::Matrix{num3}
	function GaussianCopulaNVariateProcess(rho::Matrix{num3},models::BaseProcess...) where { num3 <: Number} 
		sz=size(rho)
		@assert sz[1]==sz[2]
		@assert length(models)==sz[1]
		@assert det(rho)>=0
		return new{num3}(models,rho);
	end
	function GaussianCopulaNVariateProcess(models::BaseProcess...) 
		len_=length(models)
		return GaussianCopulaNVariateProcess(Matrix{Float64}(I, len_, len_),models...);
	end
	function GaussianCopulaNVariateProcess(model1::BaseProcess,model2::BaseProcess,rho::num3) where { num3 <: Number} 
		corr_matrix_=[1.0 rho; rho 1.0];
		@assert det(corr_matrix_)>=0
		return GaussianCopulaNVariateProcess(corr_matrix_,model1,model2);
	end
end

export GaussianCopulaNVariateProcess;

function gausscopulagen2(t::Int, Σ::Matrix{Float64} = [1. 0.5; 0.5 1.])
  z = rand(MvNormal(Σ),t)
  for i in 1:size(Σ, 2)
    d = Normal(0, sqrt(Σ[i,i]))
    @views z[i,:].= cdf.(d, z[i,:])
  end
  return collect(z');
end

function simulate!(S_total,mcProcess::GaussianCopulaNVariateProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	@assert T>0.0
	
	####Simulation
	## Simulate
	len_=length(mcProcess.models);
	for i in 1:len_
		simulate!(S_total[i],mcProcess.models[i],rfCurve,mcBaseData,T)
	end
	#S_Total::Array{Matrix{Number}}=[simulate(model_i,rfCurve,mcBaseData,T) for model_i in mcProcess.models];
	rho=mcProcess.rho
	if (rho[1,2:end]==zeros(len_-1))
		return;
	end
	for j in 1:Nstep
	
		U_joint=gausscopulagen2(Nsim,rho);
	
		for i in 1:len_
			cdf_ = sort(S_total[i][:,j+1])
			@views S_total[i][:,j+1]=Statistics.quantile(cdf_,U_joint[:,i])
		end
	end
	
	## Conclude
	return;

end
