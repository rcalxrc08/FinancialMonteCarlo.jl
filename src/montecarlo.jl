__precompile__()
module MonteCarlo

	include("pricer.jl");
	export
	    pricer,
		simulate,
		payoff



end#End Module
