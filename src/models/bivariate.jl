using EmpiricalCDFs, DatagenCopulaBased 
"""
Struct for Bidimensional Copula Processes

		gaussianCopulaBivariateProcess=GaussianCopulaBivariateProcess(σ::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		p  =	prob. of having positive jump.
		λ₊ =	positive jump size.
		λ₋ =	negative jump size.
"""
mutable struct GaussianCopulaBivariateProcess{Model_1 <: BaseProcess, Model_2 <: BaseProcess, num3 <: Number}<:BiDimensionalMonteCarloProcess
	model1::Model_1
	model2::Model_2
	rho::num3
	function GaussianCopulaBivariateProcess(σ::num,λ::num1,rho::num2) where {num <: BaseProcess, num1 <: BaseProcess, num2 <: Number}
		return new{num,num1,num2}(σ,λ,rho)
    end
end

export GaussianCopulaBivariateProcess;

function simulate(mcProcess::GaussianCopulaBivariateProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if T<=0.0
		error("Final time must be positive");
	end
	
	####Simulation
	## Simulate
	S_1=simulate(mcProcess.model1,rfCurve,mcBaseData,T);
	S_2=simulate(mcProcess.model2,rfCurve,mcBaseData,T);
	rho=mcProcess.rho
	X=[1.0 rho;rho 1.0]
	for j in 1:Nstep
		#model 1
		S_t_1_=S_1[:,j+1];
		cdf_1 = EmpiricalCDF()
		append!(cdf_1,S_t_1_)
		sort!(cdf_1)
		icdf_1 = finv(cdf_1)
		#model 2
		S_t_2_=S_2[:,j+1];
		cdf_2 = EmpiricalCDF()
		append!(cdf_2,S_t_2_)
		sort!(cdf_2)
		icdf_2 = finv(cdf_2)
		#Computation
		U_joint=gausscopulagen(Nsim,X);
		S_1[:,j+1]=icdf_1.(U_joint[:,1])
		S_2[:,j+1]=icdf_2.(U_joint[:,2])
	end
	
	## Conclude
	return (S_1,S_2);

end
