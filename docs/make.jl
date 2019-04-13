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
				"univariate.md",
				"truncate.md",
				"multivariate.md",
				"matrix.md",
				"mixture.md",
				"fit.md",
				"extends.md",
			])
deploydocs(
    repo = "https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl.git",
)