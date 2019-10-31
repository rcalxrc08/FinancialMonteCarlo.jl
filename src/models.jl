abstract type BaseProcess end

struct Underlying{num <: Number, num2 <: Number}
	S0::num
	d::num2
	function Underlying(S0::num_,d::num_2=0.0) where {num_ <: Number, num_2 <: Number}
		return new{num_,num_2}(S0,d)
	end
end
	
abstract type AbstractMonteCarloProcess <: BaseProcess end

abstract type ScalarMonteCarloProcess <: AbstractMonteCarloProcess end

abstract type VectorialMonteCarloProcess <: AbstractMonteCarloProcess end

abstract type LevyProcess<:ScalarMonteCarloProcess end

abstract type ItoProcess<:LevyProcess end

abstract type FiniteActivityProcess<:LevyProcess end

abstract type InfiniteActivityProcess<:LevyProcess end

dividend(x::mc) where {mc <: ScalarMonteCarloProcess} = mc.underlying.d;

export AbstractMonteCarloProcess
export ScalarMonteCarloProcess
export VectorialMonteCarloProcess
export ItoProcess
export LevyProcess
export FiniteActivityProcess
export InfiniteActivityProcess

using Distributions;

include("models/utils.jl")

### Ito Processes
include("models/brownian_motion.jl")
include("models/geometric_brownian_motion.jl")
include("models/black_scholes.jl")
include("models/heston.jl")
include("models/log_normal_mixture.jl")
include("models/shifted_log_normal_mixture.jl")

### Finite Activity Levy Processes
include("models/kou.jl")
include("models/merton.jl")

### Infinite Activity Levy Processes
include("models/subordinated_brownian_motion.jl")
include("models/variance_gamma.jl")
include("models/normal_inverse_gaussian.jl")


include("models/bivariate.jl")
include("models/nvariate.jl")
include("models/bivariate_log.jl")
include("models/nvariate_log.jl")

############### Display Function

import Base.Multimedia.display;
import Base.Multimedia.print;

function print(p::Union{AbstractMonteCarloProcess,AbstractPayoff})
	fldnames=collect(fieldnames(typeof(p)));
	print(typeof(p),"(");
	if(length(fldnames)>0)
		print(fldnames[1],"=",getfield(p,fldnames[1]))
		popfirst!(fldnames)
		for name in fldnames
			print(",",name,"=",getfield(p,name))
		end
	end
	print(")");
end

function display(p::Union{AbstractMonteCarloProcess,AbstractPayoff})
	fldnames=collect(fieldnames(typeof(p)));
	print(typeof(p),"(");
	print(fldnames[1],"=",getfield(p,fldnames[1]))
	popfirst!(fldnames)
	for name in fldnames
		print(",",name,"=",getfield(p,name))
	end
	println(")");
end

function get_parameters(model::XProcess) where {XProcess <: BaseProcess}
	fields_=fieldnames(XProcess)
	N1=length(fields_)
	param_=[];
	for i=1:N1
		append!(param_,getproperty(model,fields_[i]))
	end
	typecheck_=typeof(sum(param_)+prod(param_))
	param_=convert(Array{typecheck_},param_)
	
	
	return param_;
end

function set_parameters!(model::XProcess,param::Array{num}) where {XProcess <: BaseProcess , num <: Number}
	fields_=fieldnames(XProcess)
	N1=length(fields_)
	if N1!=length(param)
		error("Check number of parameters of the model")
	end
	
	for (symb,nval) in zip(fields_,param)
		setproperty!(model,symb,nval);
	end

end

function set_parameters!(model::LogNormalMixture,param::Array{num}) where {num <: Number}
	N1=div(length(param)+1,2)
	if N1*2 != length(param) + 1
		error("Check number of parameters of the model")
	end
	eta=param[1:N1]
	lam=param[(N1+1):end]
	fields_=fieldnames(LogNormalMixture)
	setproperty!(model,fields_[1],eta);
	setproperty!(model,fields_[2],lam);
end

function set_parameters!(model::ShiftedLogNormalMixture,param::Array{num}) where {num <: Number}
	N1=div(length(param),2)
	if N1*2 != length(param)
		error("Check number of parameters of the model")
	end
	eta=param[1:N1]
	lam=param[(N1+1):(2*N1-1)]
	alfa=param[end]
	fields_=fieldnames(LogNormalMixture)
	setproperty!(model,fields_[1],eta);
	setproperty!(model,fields_[2],lam);
	setproperty!(model,fields_[3],alfa);
end

"""
General Interface for Stochastic Process simulation

		S=simulate(mcProcess,spotData,mcBaseData,T)
	
Where:\n
		mcProcess          = Process to be simulated.
		spotData  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		T = Final time of the process

		S      = Matrix with path of underlying.

"""
function simulate(mcProcess::BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	error("Function used just for documentation")
end

include("models/operations.jl")