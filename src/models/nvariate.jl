"""
Struct for MultiVariate Copula Process

		gaussianCopulaNVariateProcess=GaussianCopulaNVariateProcess(models::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		models  =	the processes.
		rho  = 	correlation matrix.
"""
mutable struct GaussianCopulaNVariateProcess{ num3 <: Number, numtype <: Number} <: NDimensionalMonteCarloProcess{numtype}
	models::Tuple{Vararg{BaseProcess}}
	rho::Matrix{num3}
	function GaussianCopulaNVariateProcess(rho::Matrix{num3},models::BaseProcess...) where { num3 <: Number} 
		sz=size(rho)
		@assert sz[1]==sz[2]
		@assert length(models)==sz[1]
		@assert isposdef(rho)
		zero_typed=predict_output_type_zero_(models...)+zero(num3);
		return new{num3,typeof(zero_typed)}(models,rho);
	end
	function GaussianCopulaNVariateProcess(models::BaseProcess...) 
		len_=length(models)
		return GaussianCopulaNVariateProcess(Matrix{Float64}(I, len_, len_),models...);
	end
	function GaussianCopulaNVariateProcess(model1::BaseProcess,model2::BaseProcess,rho::num3) where { num3 <: Number} 
		corr_matrix_=[1 rho; rho 1];
		@assert isposdef(corr_matrix_)
		return GaussianCopulaNVariateProcess(corr_matrix_,model1,model2);
	end
end
support_type(z::GaussianCopulaNVariateProcess{num,num1}) where {num <: Number, num1 <: Number}=zero(num)

export GaussianCopulaNVariateProcess;

function simulate!(S_total,mcProcess::GaussianCopulaNVariateProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	@assert T>0
	
	####Simulation
	## Simulate
	len_=length(mcProcess.models);
	for i in 1:len_
		simulate!(S_total[i],mcProcess.models[i],rfCurve,mcBaseData,T)
	end
	rho=mcProcess.rho
	if (rho[1,2:end]==zeros(len_-1))
		return;
	end
	U_joint=Matrix{eltype(rho)}(undef,len_,Nsim);
	for j in 1:Nstep
	
		gausscopulagen2!(U_joint,rho,mcBaseData);
	
		for i in 1:len_
			@views tmp_=S_total[i][:,j+1]
			sort!(tmp_)
			@views tmp_.=Statistics.quantile(tmp_,U_joint[i,:];sorted=true)
		end
	end
	
	## Conclude
	return;

end
