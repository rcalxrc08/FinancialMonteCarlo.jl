# Status

_What is in flight. Layer 4: volatile, rewritten freely; GitHub issues are the source of truth for tasks, this file is
the index an agent reads first._ Updated: 2026-11-25

## Current focus

0.10.0 (2026-11-25, breaking, #72, D-014): the underlying is an argument, `pricer(model, underlying, zeroCurve, method,
payoff)`; `Underlying{S, D}` and `ZeroRate{R}` one type each (number or curve); `MarketDataSet` holds
`name → underlying` and `key → model`. Registered first in LocalRegistry for the dependent packages; the same commit goes
to General. FinancialFFT 0.8, VibratoMonteCarlo and FinancialCalibrator follow.

0.9.1 (2026-11-23): #64 (the SDE extension is triggered by StochasticDiffEq; with DifferentialEquations 8 load
StochasticDiffEq), #85 (`@non_differentiable` validations, D-013), `parameter_count` (FinancialCalibrator#36).

## Active work

| Item | Branch | Issue | Next step |
|:--|:--|:--|:--|
| Zygote on the structs | — | #80 (reopened) | constructors from vectors (#85 done: `@non_differentiable` validations, D-013) |
| New models | — | #62 | Bates, SABR, CEV (Bachelier and ShiftedProcess done) |

## Known issues

| Issue | Symptom | Workaround |
|:--|:--|:--|
| #68 | AD in `CudaMode`: only ForwardDiff works | `@test_broken` in `test/cuda/test_ad_cuda.jl` |
| #82 | arrays of payoffs/positions price an off-grid maturity at the nearest node (5-7%) | `Nstep` with the maturity on the grid |
| #74 | pathwise Greeks of barrier options biased | `@test_broken` |
| #55 | Heston biased ~10% when Feller is badly violated | — |
| #65 | Enzyme fails on a Position with legs of different payoff types | — |
| #23, #46 | Gaussian copula: path continuity, `quantile!` cost | — |

## Technical debt

- CI pins Julia 1.13.3 (Enzyme abort on 1.13.4, EnzymeAD/Enzyme.jl#4051, fixed in Enzyme 0.13.220): put '1' back once
  the Enzyme tests pass on 1.13.4.
- #66 drop the Julia 1.10 Enzyme workarounds; #67 Kou derivative in `p` with StandardMC; #83 term structure of smiles; #28 European payoffs without the whole path
  matrix.
- Untracked local files: `lit.jl`, `lit2.jl`, `path.jl`, `examples/`, `benchmark/bench_black_ql.jl`, `benchmark/wip/`,
  `src/02_models/models/ivm/`, `src/03_options/options/payoffs/basket/n_average_forward.jl`, `test/wip/` additions; an
  uncommitted change in `test/models/diffeq.jl`.

## Working notes

Notes of an active initiative live next to this file as `state/YYYY-MM-DD-<topic>.md`. When the initiative ends, move
what stays true into `references/` or `references/decisions.md` and delete the note.
