__precompile__()

module MonteCarlo
	DIFFEQ_MONTECARLO_ACTIVE_FLAG=true
	if(!((VERSION.major==0)&&(VERSION.minor<=6)))
		using Random, LinearAlgebra
	end
	include("utils.jl")
	include("payoffs.jl")
	include("models.jl")
	include("metrics.jl");
	export
	    pricer,
	    variance,
		confinter,
		simulate,
		payoff,
		DIFFEQ_MONTECARLO_ACTIVE_FLAG



end#End Module
