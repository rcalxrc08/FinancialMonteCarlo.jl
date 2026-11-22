# Documentation

The Documenter + Literate site (`Documentation.yml`, GitHub Pages from Actions, master at the root, the latest release
in `stable/`).

## Inputs

| Load | Why |
|:--|:--|
| `make.jl` | pages in order (tutorials, manual `Manual NN`, extending, API), Literate rendering, no `deploydocs` |
| `figures.jl` | figures simulated at build time (light and dark SVG), `hero_figure` |
| `literate/getting_started.jl`, `literate/exotics_and_greeks.jl`, `literate/portfolios.jl` | the tutorials |
| `references/docs_design.md` | the shared design, the manual structure, the bestiaries |
| `references/workflows.md` | the deployment of master and `stable/` |
| `references/public_api.md` | what `api.md` lists |

## Process

- Keep the manual, `api.md`, `extending.md`, the tutorials and the docstrings in step with the code in the same commit.
- A new model gets a card in `src/manual/bestiary.md`; a new payoff a row in `src/manual/payoff_bestiary.md`.
- Measured numbers say on which machine they were taken.

## Outputs

A site that builds with no warnings.

## Verify (Kaimon)

| Check | Expect |
|:--|:--|
| Docs build in a separate process: `julia --project=docs docs/make.jl` | builds, no warnings |
