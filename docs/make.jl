# Build the documentation site:
#     julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#     julia --project=docs docs/make.jl
# The HTML ends up in docs/build; on GitHub Actions (.github/workflows/Documentation.yml) it is uploaded as the Pages
# artifact and deployed with GitHub Pages, at the site root (master) and in stable/ (the latest release).
using Documenter, Literate, FinancialMonteCarlo

# The tutorials draw plots with GR, which needs no display in this mode (CI runners have none).
ENV["GKSwstype"] = "100"

# Tutorials are Literate.jl scripts in docs/literate, rendered to markdown (and executed) at build time.
const LITERATE_DIR = joinpath(@__DIR__, "literate")
const TUTORIALS_DIR = joinpath(@__DIR__, "src", "tutorials")
const TUTORIALS = ["getting_started.jl", "exotics_and_greeks.jl", "portfolios.jl"]
for tutorial in TUTORIALS
    Literate.markdown(joinpath(LITERATE_DIR, tutorial), TUTORIALS_DIR; documenter = true)
end

# The figures of the home page are computed with the package, in a light and a dark variant (see figures.jl and
# assets/figures.css), in a module of their own.
module DocumentationFigures
include(joinpath(@__DIR__, "figures.jl"))
end

makedocs(
    sitename = "FinancialMonteCarlo.jl",
    authors = "rcalxrc08",
    modules = [FinancialMonteCarlo],
    repo = Remotes.GitHub("rcalxrc08", "FinancialMonteCarlo.jl"),
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/favicon.ico", "assets/figures.css", "assets/style.css", "assets/site.css", "assets/home.css", "assets/dark.js"], canonical = "https://rcalxrc08.github.io/FinancialMonteCarlo.jl/"),
    pages = ["Home" => "index.md", "Tutorials" => ["Getting started" => "tutorials/getting_started.md", "Exotic options and greeks" => "tutorials/exotics_and_greeks.md", "Portfolios and multivariate models" => "tutorials/portfolios.md"], "Manual" => ["Overview" => "manual/overview.md", "Market data" => "manual/market_data.md", "Monte Carlo configuration" => "manual/configuration.md", "Models" => "manual/models.md", "Payoffs" => "manual/payoffs.md", "Metrics" => "manual/metrics.md", "Portfolios and multivariate" => "manual/portfolios.md", "Automatic differentiation" => "manual/autodiff.md", "DifferentialEquations.jl" => "manual/diffeq.md", "Performance" => "manual/performance.md", "Extending the package" => "extending.md"], "Reference" => ["API" => "api.md"]],
    warnonly = [:missing_docs]
)
