__precompile__()

module MonteCarlo
	DIFFEQ_MONTECARLO_ACTIVE_FLAG=true
	if(!((VERSION.major==0)&&(VERSION.minor<=6)))
		using Random, LinearAlgebra
	end
	using Requires
	function __init__()
		@require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("cuda_dependencies.jl")
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
