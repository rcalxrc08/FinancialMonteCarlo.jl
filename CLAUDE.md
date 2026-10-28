# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FinancialMonteCarlo.jl: a Julia package for Monte Carlo pricing and risk of equity derivatives (single-name and multivariate via Gaussian copula). The canonical remote is GitLab (`gitlab.com/rcalxrc08/FinancialMonteCarlo.jl`, docs deployed via GitLab Pages); GitHub Actions also run CI (Julia 1.10 and latest 1.x), CompatHelper and TagBot.

## Commands

```bash
# Run the full test suite (add coverage = true for .cov files)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run only some groups (names from TEST_GROUPS in test/runtests.jl, case-insensitive)
julia --project=. -e 'using Pkg; Pkg.test(test_args = ["models", "ad"])'
julia --project=test test/runtests.jl models

# Run a single test file (test/Project.toml points at the package through [sources])
julia --project=test -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=test test/models/black_scholes.jl

# Benchmarks
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("benchmark/benchmarks.jl")'

# Docs (Documenter + Literate; docs/literate/*.jl are rendered into docs/src/tutorials at build time)
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl

# Formatting (JuliaFormatter, config in .JuliaFormatter.toml — margin=1200, so long one-line signatures are expected)
julia -e 'using JuliaFormatter; format(".")'
```

### Test harness specifics

- [test/runtests.jl](test/runtests.jl) runs an **explicit** list of files grouped in `TEST_GROUPS` (Structs, Models, Methods, Metrics, AD, Portfolio); a new test file must be added there. Each file is included into its own module, so files cannot share globals.
- Every test file starts with `include("../testutils.jl")` (or `"testutils.jl"` at top level) so it also runs standalone. [test/testutils.jl](test/testutils.jl) holds the shared market data (`S0, K, r, T, d, D, sigma`), the configurations `mc` / `mc_anti`, reference prices (`BS_CALL_REF`, `BS_PUT_REF`), `standard_payoffs()` and `check_prices(model, rfCurve, mc, payoffs, refs; toll)` (one testset per payoff: compares with `refs` when present, otherwise checks the price is finite). `capture_stdout` tests the package's `display`/`print` overloads, which write to stdout directly.
- Tests compare seeded MC prices with reference values within an absolute tolerance, so a change in RNG consumption can shift prices without being a bug.
- `test/cuda/` (needs a GPU) and `test/wip/` (scratch, outdated API) are not run.
- [benchmark/benchmarks.jl](benchmark/benchmarks.jl) still includes every file in `benchmark/` except a blacklist.

## Architecture

Entry point: [src/FinancialMonteCarlo.jl](src/FinancialMonteCarlo.jl). Include order matters (`utils` → `models` → `options` → `io_utils` → `operations` → `metrics` → `output_type` → `multi_threading`). Each type file does its own `export`.

### Core pipeline

Every metric reduces to: `pricer(process, zeroCurve, mcConfig, payoff)` → `set_seed!` → `simulate(process, zeroCurve, mcConfig, maturity(payoff))` returns an `Nsim × (Nstep+1)` path matrix → `payoff(S, payoff, zeroCurve, mcConfig, T)` returns discounted per-path payoffs → `mean`. See [src/metrics/pricer.jl](src/metrics/pricer.jl). `pricer` is overloaded for: a single payoff, an array of payoffs, a `Position`, an `AbstractSmile`, a vectorial process with `Array{Position}`, and a `MarketDataSet` with a `Portfolio`. `delta`, `rho`, `variance`, `confinter`, `distribution` build on top (delta/rho are bump-and-reprice using `Accessors.@set` on the immutable model structs).

### Models ([src/models/](src/models/))

Type hierarchy in [base_models.jl](src/models/base_models.jl):
- `AbstractMonteCarloEngine` — raw engines with no notion of rates/dividends (`BrownianMotion`, `GeometricBrownianMotion`, `SubordinatedBrownianMotion`); `simulate!(X, engine, mcConfig, T)`.
- `AbstractMonteCarloProcess` — financial models holding an `underlying::AbstractUnderlying` (S0, dividend); `simulate!(X, process, zeroCurve, mcConfig, T)`. Subtypes: `ItoProcess` (BlackScholes, Heston, LogNormalMixture…), `FiniteActivityProcess` (Kou, Merton), `InfiniteActivityProcess` (VarianceGamma, NIG), and `VectorialMonteCarloProcess`/`NDimensionalMonteCarloProcess` (copula-based `GaussianCopulaNVariateProcess` in `nvariate*.jl`).
- Models are typically thin: they compute a risk-neutral drift and delegate to an engine, then exponentiate. E.g. [fa_levy.jl](src/models/fa_levy.jl) implements `simulate!` generically for all finite-activity Lévy models; concrete models (Kou, Merton) only provide `compute_jump_size` and `compute_drift`.

Variants are selected **by multiple dispatch, not by separate types**:
- `*_aug.jl` files: methods for non-constant curves (`ZeroRateCurve` / `CurveType` drift, e.g. `BrownianMotionVec`).
- `*_sobol.jl` files: methods dispatching on `SerialSobolMonteCarloConfig`.
- Antithetic: methods on `SerialAntitheticMonteCarloConfig` (requires even `Nsim`).
These type aliases are defined in [src/utils/mc_config.jl](src/utils/mc_config.jl).

### Configuration ([src/utils/](src/utils/))

`MonteCarloConfiguration(Nsim, Nstep, method, parallelMode)`:
- `method <: AbstractMonteCarloMethod`: `StandardMC`, `AntitheticMC`, `SobolMode`, `PrescribedMC`, or `ControlVariates` (defined in [payoff_base.jl](src/options/payoff_base.jl), handled in [pricer_cv.jl](src/metrics/pricer_cv.jl)).
- `parallelMode <: BaseMode`: `SerialMode(seed, rng)` (models draw from `mcConfig.parallelMode.rng`), `MultiThreading` ([src/multi_threading.jl](src/multi_threading.jl) — splits Nsim into seeded batches and averages; pricer methods generated via `@eval`), plus optional CUDA/Distributed modes.
- Rates: `ZeroRate(r)` scalar or `ZeroRateCurve`/`ImpliedZeroRate` built on `CurveType` ([curves.jl](src/utils/curves.jl)); code uses `integral(r, T)` / `incremental_integral` so the same path works for scalar and curve rates.

### Payoffs ([src/options/](src/options/))

Hierarchy in [payoff_base.jl](src/options/payoff_base.jl): `SingleNamePayoff` → `EuropeanPayoff`, `AmericanPayoff`, `BermudanPayoff`, `PathDependentPayoff` (→ `BarrierPayoff`, `AsianPayoff`); separately `BasketPayoff`. Generic engines in [src/options/engines/](src/options/engines/) implement `payoff(S, ...)` per family (American/Bermudan use Longstaff–Schwartz); a concrete payoff usually only defines its struct and a scalar `payout(ST, payoff)` (European) or path-based equivalent. Payoffs must expose `.T` (used by `maturity`). Smiles ([src/options/smiles/](src/options/smiles/)) are strike-vectors expanded via `smile_to_options`.

### Contract algebra / portfolios ([src/operations/](src/operations/))

- `Position = Dict{AbstractPayoff, Number}`: payoffs combine with `+`, `*`, `/` into weighted positions.
- `MarketDataSet = Dict{String, AbstractMonteCarloProcess}` built with `"KEY" → model` (the `→` operator is exported); multivariate keys use `"IDX1_IDX2"`.
- `Portfolio` maps underlying keys to positions; `pricer(::MarketDataSet, ..., ::Portfolio)` simulates each underlying once up to the max maturity.

### Automatic differentiation conventions

AD (ForwardDiff, ReverseDiff, Mooncake, DualNumbers, TaylorSeries, HyperDualNumbers) is a first-class concern, which drives several pervasive patterns — follow them in new code:
- Structs are fully parametric on each numeric field and carry a promoted `numtype` parameter (e.g. `KouProcess{num, num1, …, numtype}`), computed in the inner constructor via `promote_type` or `zero(num1) + zero(num2)`.
- Output element types are inferred with `predict_output_type_zero(...)` ([src/output_type.jl](src/output_type.jl)) to allocate path/price arrays of the right (possibly Dual) type — never hardcode `Float64` buffers.
- Assertions/validation and non-differentiable bookkeeping are wrapped in `ChainRulesCore.@ignore_derivatives` (needed for Zygote/ReverseDiff compatibility).
- Mooncake needs a pure-Julia RNG (`SerialMode(seed, Xoshiro())`): the default `MersenneTwister` calls into libdSFMT, which Mooncake cannot differentiate. Mooncake's first rule build is slow (~2 min in `test/ad/mooncake.jl`).

### Optional dependencies

[src/deps/](src/deps/) files are loaded lazily through `Requires.@require` in `__init__`: CUDA (`*.cujl` files, GPU models in [src/models/cuda/](src/models/cuda/)), DifferentialEquations (price directly from an SDE via [diff_eq_monte_carlo.jl](src/models/diff_eq_monte_carlo.jl)), DualNumbers, TaylorSeries, VectorizedRNG, Distributed. Code there only exists when the user has loaded the corresponding package.

## Extending (see docs/src/extending.md for runnable examples)

- New model: struct subtyping the appropriate abstract type (engine vs. process matters), holding its parameters (+ `underlying` for processes), and a `simulate!` method; copy the signature from an existing model.
- New payoff: struct subtyping the right payoff family with its parameters (`T`, `K`, barriers…) and a `payoff`/`payout` method; all metrics then work automatically.
- `get_parameters` / `set_parameters!` exist on models and payoffs and are tested — keep them consistent when adding fields.

## Caveats (from README)

- First- and second-order AD "works" for every payoff, including ones where it is meaningless (binary, barrier); second-order results are mostly useless.
- No multicurrency support.
