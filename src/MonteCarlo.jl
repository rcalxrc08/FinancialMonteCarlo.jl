__precompile__()

module MonteCarlo
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
		payoff



end#End Module
