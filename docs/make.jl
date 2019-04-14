using Documenter, FinancialMonteCarlo

makedocs(
		format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
		sitename="FinancialMonteCarlo.jl",
		modules = [FinancialMonteCarlo],
		pages = [
				"index.md",
				"starting.md",
				"types.md",
				"stochproc.md",
				"parallel_vr.md",
				"multivariate.md",
				"payoffs.md",
				"mixture.md",
				"fit.md",
				"extends.md",
			])
deploydocs(
    repo = "https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl.git",
)