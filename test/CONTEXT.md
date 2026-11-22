# Tests and CI

The groups, the shared data, what CI runs.

## Inputs

| Load | Why |
|:--|:--|
| `runtests.jl` | `TEST_GROUPS` (Structs, Models, Methods, Metrics, AD, Portfolio); each file in its own module, logged |
| `testutils.jl` | market data, `mc`/`mc_anti`/`mc_small`/`mc_small_anti`, `BS_CALL_REF`/`BS_PUT_REF`, `standard_payoffs`, `check_prices`, `capture_stdout` |
| `structs.jl`, `engines.jl`, `aqua.jl`, `context.jl` | Structs group: constructors and checks, Longstaff–Schwartz against a BigFloat fit, Aqua, the context layers |
| `models/black_scholes.jl`, `models/barrier_closed_form.jl`, `models/bachelier.jl`, `models/heston.jl`, `models/kou.jl`, `models/merton.jl`, `models/variance_gamma.jl`, `models/normal_inverse_gaussian.jl`, `models/lognormal_mixture.jl`, `models/shifted.jl`, `models/parameters.jl`, `models/diffeq.jl` | Models group |
| `methods/sobol.jl`, `methods/control_variates.jl`, `methods/parallel.jl`, `methods/multiprocess.jl` | Methods group |
| `metrics/statistics.jl`, `metrics/greeks.jl` | Metrics group (Float32 in, Float32 out) |
| `ad/dual_numbers.jl`, `ad/forward_diff.jl`, `ad/reverse_diff.jl`, `ad/mooncake.jl`, `ad/enzyme.jl`, `ad/taylor_series.jl`, `ad/chain_rules.jl` | AD group |
| `portfolio/algebra.jl`, `portfolio/multivariate.jl` | Portfolio group |
| `cuda/test_ad_cuda.jl`, `cuda/test_black_cuda.jl`, `cuda/test_device_paths_cuda.jl`, `cuda/test_fa_levy_cuda.jl`, `cuda/test_new_models_cuda.jl`, `cuda/test_smile_cuda.jl` | local GPU tests, not run by `runtests.jl` |
| `references/testing.md` | reference prices and their sources, durations, what is not run |
| `references/workflows.md` | commands, CI matrix |

## Process

- A new test file is added to `TEST_GROUPS` and starts with `include("../testutils.jl")`.
- Reference prices from FinancialToolbox closed forms; seeded comparisons within an absolute tolerance (a change in RNG
  consumption can shift prices without being a bug).
- `test/wip/` is scratch with an outdated API, not run.

## Outputs

A green `Pkg.test()` (~20 minutes locally); `context.jl` green whenever a file is added, moved or renamed.

## Verify (Kaimon)

| Check | Expect |
|:--|:--|
| `run_tests()` summary | pass; reports the Julia of `kaimon.toml` |
| `include("test/context.jl")` | pass |
