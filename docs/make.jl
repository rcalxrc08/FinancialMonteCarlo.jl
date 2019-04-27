using Documenter, FinancialMonteCarlo

makedocs(
		format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
		assets = ["assets/favicon.ico"],
		assets = ["assets/logo.png"]
    ),
		sitename="FinancialMonteCarlo.jl",
		modules = [FinancialMonteCarlo],
		pages = [
				"index.md",
				"starting.md",
				"types.md",
				"stochproc.md",
				"parallel_vr.md",
				"payoffs.md",
				"intdiffeq.md",
				"extends.md",
			])
get(ENV, "CI", nothing) == "true" ? deploydocs(
    repo = "https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl.git",
) : nothing