# Testing

_The harness, the shared data, reference values, durations._ (Layer 3 reference.)

- `test/runtests.jl` runs an **explicit** list of files grouped in `TEST_GROUPS` (Structs, Models, Methods, Metrics, AD,
  Portfolio); a new test file must be added there. Each file is included into its own module, so files cannot share
  globals. The testsets report only at the end, so `run_test_file` logs `Running <file>` and `Finished <file> in N s`
  (flushed): a hung file shows up by name in a local log or on CI.
- A full local run takes about 20 minutes; run the groups one by one to see a summary per group. Most of the duration is
  compilation (`test/ad/enzyme.jl` alone is about 2.5 minutes, `test/ad/mooncake.jl` 1.5), so a new combination of
  model, payoff and numeric type costs more than its `Nsim`.
- `test/aqua.jl` (Structs group) runs Aqua.jl: method ambiguities, type piracy, unbound type parameters, stale
  dependencies and compat bounds must stay clean. `test/context.jl` (Structs group) checks the context layers. For type
  stability, JET's `report_opt` on `pricer`/`variance` finds only the intended runtime dispatch at the Sobol function
  barrier (`on_sobol_paths`).
- Every test file starts with `include("../testutils.jl")` (or `"testutils.jl"` at top level) so it also runs
  standalone. `test/testutils.jl` holds the shared market data (`S0, K, r, T, d, D, sigma`), the configurations `mc` /
  `mc_anti` (for comparisons with reference prices) and `mc_small` / `mc_small_anti` (for what does not depend on the
  Monte Carlo error: shapes, finite prices, identities on the same seed), reference prices (`BS_CALL_REF`,
  `BS_PUT_REF`: the forward, European and binary European entries are the closed forms of FinancialToolbox,
  `blsprice`/`blsbin`, the American, barrier and Asian ones hard-coded estimates), `standard_payoffs()` and
  `check_prices(model, rfCurve, mc, payoffs, refs; toll)` (one testset per payoff: compares with `refs` when present,
  otherwise checks the price is finite). `capture_stdout` tests the package's `display`/`print` overloads.
- `test/models/barrier_closed_form.jl` checks the four single barriers against `blsdownout`/`blsdownin`/`blsupout`/
  `blsupin` of FinancialToolbox ≥ 0.8.1 with the Broadie-Glasserman-Kou shift of the barrier for discrete monitoring,
  and records the bias of their pathwise delta as `@test_broken` (#74).
- Tests compare seeded MC prices with reference values within an absolute tolerance, so a change in RNG consumption can
  shift prices without being a bug.
- `test/cuda/` (needs a GPU) and `test/wip/` (scratch, outdated API) are not run. `test/cuda/test_ad_cuda.jl` checks the
  AD backends in `CudaMode` against finite differences (needs Enzyme and Mooncake in the environment, about 3 hours);
  the upstream causes are in the comment above `BROKEN` in that file (#68, `autodiff.md`).
- Closed forms in tests and docs come from FinancialToolbox (`blsprice`, `blsbin`, `blsimpv`), never written by hand; it
  is a dependency of `docs/` and `test/`. The Bachelier closed form waits for rcalxrc08/FinancialToolbox.jl#84.
