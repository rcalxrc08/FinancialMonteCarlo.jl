# Build the documentation site:
#     julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#     julia --project=docs docs/make.jl
# The HTML ends up in docs/build (published as GitLab Pages by .gitlab-ci.yml).
using Documenter, Literate, FinancialMonteCarlo

# Tutorials are Literate.jl scripts in docs/literate, rendered to markdown (and executed) at build time.
const LITERATE_DIR = joinpath(@__DIR__, "literate")
const TUTORIALS_DIR = joinpath(@__DIR__, "src", "tutorials")
const TUTORIALS = ["getting_started.jl", "exotics_and_greeks.jl", "portfolios.jl"]
for tutorial in TUTORIALS
    Literate.markdown(joinpath(LITERATE_DIR, tutorial), TUTORIALS_DIR; documenter = true)
end

makedocs(
    sitename = "FinancialMonteCarlo.jl",
    authors = "rcalxrc08",
    modules = [FinancialMonteCarlo],
    repo = Remotes.GitLab("rcalxrc08", "FinancialMonteCarlo.jl"),
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/favicon.ico"], canonical = "https://rcalxrc08.gitlab.io/FinancialMonteCarlo.jl/"),
    pages = [
        "Home" => "index.md",
        "Tutorials" => ["Getting started" => "tutorials/getting_started.md", "Exotic options and greeks" => "tutorials/exotics_and_greeks.md", "Portfolios and multivariate models" => "tutorials/portfolios.md"],
        "Manual" => [
            "Overview" => "manual/overview.md",
            "Market data" => "manual/market_data.md",
            "Monte Carlo configuration" => "manual/configuration.md",
            "Models" => "manual/models.md",
            "Payoffs" => "manual/payoffs.md",
            "Metrics" => "manual/metrics.md",
            "Portfolios and multivariate" => "manual/portfolios.md",
            "Automatic differentiation" => "manual/autodiff.md",
            "DifferentialEquations.jl" => "manual/diffeq.md",
        ],
        "Extending the package" => "extending.md",
        "API reference" => "api.md",
    ],
    # Docstrings of the optional (Requires-based) extensions are not loaded here.
    warnonly = [:missing_docs],
)
