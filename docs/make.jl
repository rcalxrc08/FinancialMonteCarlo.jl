using Documenter, FinancialMonteCarlo, Literate
getting_started = joinpath(@__DIR__, "..", "examples", "getting_started.jl")
intdiffeq = joinpath(@__DIR__, "..", "examples", "intdiffeq.jl")
OUTPUT = joinpath(@__DIR__, "src")
Literate.markdown(getting_started, OUTPUT; documenter = true)
Literate.markdown(intdiffeq, OUTPUT; documenter = true)
makedocs(format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/favicon.ico"]), sitename = "FinancialMonteCarlo.jl", modules = [FinancialMonteCarlo], pages = ["index.md", "getting_started.md", "types.md", "stochproc.md", "parallel_vr.md", "payoffs.md", "metrics.md", "multivariate.md", "intdiffeq.md", "extends.md"])
get(ENV, "CI", nothing) == "true" ? deploydocs(repo = "https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl.git") : nothing