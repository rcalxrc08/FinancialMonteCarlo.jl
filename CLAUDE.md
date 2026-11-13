# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FinancialMonteCarlo.jl: a Julia package for Monte Carlo pricing and risk of equity derivatives (single-name and multivariate via Gaussian copula). The canonical remote is GitLab (`gitlab.com/rcalxrc08/FinancialMonteCarlo.jl`); GitHub Actions run CI (Julia 1.10 and latest 1.x), the benchmarks, the documentation (built by [Documentation.yml](.github/workflows/Documentation.yml) and published with GitHub Pages from Actions, with no `gh-pages` branch: master at the site root, the latest release in `stable/`; GitHub Pages does not serve a deployment made from a tag, so a tag run starts a run on master, which rebuilds `stable/` from the new tag), CompatHelper and TagBot.

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

# Benchmarks (benchmark/Project.toml points at the package through [sources]); prints a markdown table at the end
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=benchmark -t auto benchmark/run.jl                 # all groups
BENCHMARK_SECONDS=0.3 julia --project=benchmark -t auto benchmark/run.jl payoffs greeks   # some groups, quickly

# Docs (Documenter + Literate; docs/literate/*.jl are rendered into docs/src/tutorials at build time)
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl

# Formatting (JuliaFormatter, config in .JuliaFormatter.toml — margin=1200, so long one-line signatures are expected)
julia -e 'using JuliaFormatter; format(".")'
```

### Test harness specifics

- [test/runtests.jl](test/runtests.jl) runs an **explicit** list of files grouped in `TEST_GROUPS` (Structs, Models, Methods, Metrics, AD, Portfolio); a new test file must be added there. Each file is included into its own module, so files cannot share globals. `test/aqua.jl` (Structs group) runs Aqua.jl: method ambiguities, type piracy, unbound type parameters, stale dependencies and compat bounds must stay clean. For type stability, JET's `report_opt` on `pricer`/`variance` finds only the intended runtime dispatch at the Sobol function barrier (`on_sobol_paths`).
- Every test file starts with `include("../testutils.jl")` (or `"testutils.jl"` at top level) so it also runs standalone. [test/testutils.jl](test/testutils.jl) holds the shared market data (`S0, K, r, T, d, D, sigma`), the configurations `mc` / `mc_anti` (for comparisons with reference prices) and `mc_small` / `mc_small_anti` (for what does not depend on the Monte Carlo error: shapes, finite prices, identities on the same seed), reference prices (`BS_CALL_REF`, `BS_PUT_REF`), `standard_payoffs()` and `check_prices(model, rfCurve, mc, payoffs, refs; toll)` (one testset per payoff: compares with `refs` when present, otherwise checks the price is finite). `capture_stdout` tests the package's `display`/`print` overloads, which write to stdout directly.
- Tests compare seeded MC prices with reference values within an absolute tolerance, so a change in RNG consumption can shift prices without being a bug.
- Most of the suite's duration is compilation (`test/ad/enzyme.jl` alone is about 2.5 minutes, `test/ad/mooncake.jl` 1.5), so a new combination of model, payoff and numeric type costs more than its `Nsim`.
- `test/cuda/` (needs a GPU) and `test/wip/` (scratch, outdated API) are not run.
- [benchmark/benchmarks.jl](benchmark/benchmarks.jl) defines `SUITE`, a BenchmarkTools group per area (`models`, `methods`, `payoffs`, `markets`, `greeks`; 10 000 paths × 50 steps); [benchmark/run.jl](benchmark/run.jl) tunes and runs the groups given as arguments (all by default, `BENCHMARK_SECONDS` per benchmark, 2 s by default) and prints a markdown table (median, minimum, memory, allocations), also written to the job summary on GitHub. [Benchmarks.yml](.github/workflows/Benchmarks.yml) runs it on every push and pull request to master and by hand (`gh workflow run Benchmarks.yml -f groups="payoffs"`); shared runners vary by 10-20%, so compare commits within one machine. The whole suite takes about 4 minutes there (first run, 5eb3e1c: AMD EPYC 7763, 4 threads, European on Black-Scholes 4.1 ms). A new feature with a cost of its own gets a benchmark there.

## Architecture

Entry point: [src/FinancialMonteCarlo.jl](src/FinancialMonteCarlo.jl). Include order matters (`utils` → `models` → `options` → `io_utils` → `operations` → `metrics` → `output_type` → `batches` → `multi_threading`). Each type file does its own `export`.

### Core pipeline

Every metric reduces to: `pricer(process, zeroCurve, mcConfig, payoff)` → `set_seed!` → `simulate(process, zeroCurve, mcConfig, maturity(payoff))` returns an `Nsim × (Nstep+1)` path matrix → `payoff(S, payoff, zeroCurve, mcConfig, T)` returns discounted per-path payoffs → `mean`. See [src/metrics/pricer.jl](src/metrics/pricer.jl). `pricer` is overloaded for: a single payoff, an array of payoffs, a `Position`, an `AbstractSmile`, a vectorial process with `Array{Position}`, and a `MarketDataSet` with a `Portfolio`. `delta`, `rho`, `variance`, `confinter`, `distribution` build on top (delta/rho are bump-and-reprice using `Accessors.@set` on the immutable model structs).

### Models ([src/models/](src/models/))

Type hierarchy in [base_models.jl](src/models/base_models.jl):
- `AbstractMonteCarloEngine` — raw engines with no notion of rates/dividends (`BrownianMotion`, `GeometricBrownianMotion`, `SubordinatedBrownianMotion`); `simulate!(X, engine, mcConfig, T)`.
- `AbstractMonteCarloProcess` — financial models holding an `underlying::AbstractUnderlying` (S0, dividend); `simulate!(X, process, zeroCurve, mcConfig, T)`. Subtypes: `ItoProcess` (BlackScholes, Heston, LogNormalMixture…), `FiniteActivityProcess` (Kou, Merton), `InfiniteActivityProcess` (VarianceGamma, NIG), and `VectorialMonteCarloProcess`/`NDimensionalMonteCarloProcess` (copula-based `GaussianCopulaNVariateProcess` in `nvariate*.jl`).
- Models are typically thin: they compute a risk-neutral drift and delegate to an engine, then exponentiate. E.g. [fa_levy.jl](src/models/fa_levy.jl) implements `simulate!` generically for all finite-activity Lévy models; concrete models (Kou, Merton) only provide `compute_jump_size` and `compute_drift`. Likewise [ia_levy.jl](src/models/ia_levy.jl) for infinite-activity Lévy models: concrete models (VarianceGamma, NIG) only provide `compute_drift` and `compute_subordinator_law`, each in two forms: one reading the model's fields, delegating to one taking the parameters explicitly (`compute_drift(model, σ, θ, κ)`, `compute_subordinator_law(model, κ, dt)`), which can be called with other values and numeric types, e.g. to build a subordinator without Duals.

A model has **one `simulate!` body** for antithetic variates, constant and curve rates, and CPU and GPU, built on the hooks in [src/models/paths.jl](src/models/paths.jl), so do not add a copy per variant:
- `step_drifts(μ, dt, Nstep)` / `drift_of_step(drifts, j)`: a constant rate stays one `Number` and only a curve is a vector, so flat simulations allocate nothing for the drift.
- `path_sides(X, mcConfig)`: all the paths or, with `AntitheticMC`, the two halves with the sign of their normal draws.
- `draw_buffer`, `draw_normals!`, `draw_subordinator!`, `get_matrix_type`: buffers and draws on the device of the mode. The arrays come from `array_type(mode)`; a device adds its method of `array_type` and, for the arrays of its own type (`Z::CuArray`), of `draw_normals!` and `draw_subordinator!(dt_s, law, mcConfig)` ([ext/FinancialMonteCarloCUDAExt/](ext/FinancialMonteCarloCUDAExt/)). Models draw on the host from `host_rng(mcConfig)`.
- What CUDA keeps off the host: the `Gamma` and `InverseGaussian` subordinators are drawn by a kernel (`device_draw`, Marsaglia-Tsang and Michael-Schucany-Haas; any other law is sampled on the host by Distributions and copied), whose generator is seeded at every launch from `host_rng`, because `CUDA.seed!` and `Random.seed!` do not reach the generator of the kernels. The jumps of Kou and Merton are drawn on the host into a `JumpList` (the samplings write through `add_jump!(sink, i, jump_idx, size)`, the sink being the matrix or the list) and added by one kernel with a thread per element (contiguous writes; a thread per path writing along its row was 4 times slower). Barrier and Asian payoffs are row reductions (`payoff_by_rows`, from `barrier_payout(ST, lo, hi, payoff)` / `asian_payout(ST, A, payoff)`; `is_alive` takes the path or its `PathRange`), so a new barrier option needs nothing and a new Asian one its `asian_payout` and a place in the `Union` of the extension.

Sobol is the one separate body (`*_sobol.jl`, dispatching on `SobolConfig = MonteCarloConfiguration{<:SobolMC}`), because it draws path by path (the first `2^m - 1` points are skipped by `sobol_skip!`, which sets the state of the sequence in closed form: never call `Sobol.skip`, which draws every skipped point); it fills each path inside `on_sobol_paths`, a function barrier, because the type of a `SobolSeq` holds its dimension, known only at run time. Kou and Merton (`FiniteActivityProcess`) draw their jumps from uniforms by inversion with `AntitheticMC` and `SobolMC` ([fa_levy.jl](src/models/fa_levy.jl)): the count from the Poisson table `JumpUniforms` (`jump_count`), then a time and a size (`jump_quantile`, one method per model) per jump; antithetic paths take `1 - u`, and a Sobol point has `Nstep + 1 + 2M` coordinates (`on_sobol_paths(...; normals = Nstep)` leaves the jump coordinates uniform), beyond `M` jumps the rng. VibratoMonteCarlo reuses these helpers. The aliases `SobolConfig` and `BatchedConfig` are defined in [src/utils/mc_config.jl](src/utils/mc_config.jl).

### Configuration ([src/utils/](src/utils/))

Three separate concepts, in [mc_mode.jl](src/utils/mc_mode.jl) and [mc_config.jl](src/utils/mc_config.jl):
- **Method** (`AbstractMethod`): what `pricer(model, zeroCurve, method, payoff)` dispatches on. `MonteCarloConfiguration <: AbstractMonteCarloConfiguration <: AbstractMethod` here; FinancialFFT, VibratoMonteCarlo and FinancialCalibrator subtype it too. **Every method carries its execution mode**: `execution_mode(method)` reads it, `with_mode(method, mode)` returns a copy in another mode (generic `@set method.mode`; `MonteCarloConfiguration` rebuilds its random state). Never add a trailing `mode` argument to a metric.
- **Mode** (`AbstractMode`): where and how a method runs, with no seed and no generator, so methods without randomness use the same objects. `SerialMode()`, `MultiThreading(nbatches)`, `MultiProcess(nbatches)` (batched modes, `AbstractBatchedMode`), `CudaMode()` (`AbstractCudaMode`, defined in `src/` so other packages can dispatch on it without CUDA).
- **Sampling** (`AbstractSampling`): `StandardMC`, `AntitheticMC`, `SobolMC`, or `ControlVariates` (defined in [payoff_base.jl](src/options/payoff_base.jl), handled in [pricer_cv.jl](src/metrics/pricer_cv.jl)).

`MonteCarloConfiguration(Nsim, Nstep; sampling = StandardMC(), mode = SerialMode(), seed = 0, rng = MersenneTwister(), reseed = true)` is the only public constructor (fields `Nsim, Nstep, sampling, mode, rng, reseed`; the default 6-argument constructor is internal, used by `batch_configs` and `@set`). It runs every compatibility check once, in `check_configuration` (`check_sampling`, `check_sampling_mode`, `check_reseed`, `check_mode`), so an unsupported combination fails at construction; add new restrictions there as a dispatch method, not in the metrics. Its `rng` field is the random state built by `random_state(mode, seed, rng)`: one `SeededRNG(seed, rng)`, or for a batched mode a `Vector{SeededRNG}` (batch `i` gets a copy of `rng` with seed `seed + i`, or `seed[i]` when `seed` is a vector). `set_seed!` reseeds them (if `reseed`) and calls `seed_device!(mode, seed)` (CURAND in the CUDA extension). Batched modes split into serial batch configurations ([src/batches.jl](src/batches.jl) holds all they share, and a mode only adds a `run_on_batches` method, whose results `batch_results` converts to the type inference gives (`fetch` and `pmap` return `Any`, which would make every batched metric type unstable): `pricer` averages the batches, `variance`/`confinter` pool their statistics, `simulate` joins their paths), so results do not depend on threads or workers. The 0.5 names (`SobolMode`, `BaseMode`, `AbstractMonteCarloMethod`, positional constructors) are deprecated in [src/deprecated.jl](src/deprecated.jl).

**Data placement**: the mode decides where the computation runs and where results live. Inputs holding arrays are moved to the mode once when a metric starts, with `to_mode(x, mode)` (`Adapt.adapt(mode_adaptor(mode), x)`; `HostAdaptor` for the CPU modes, `CuArray` for CUDA). Smiles define `Adapt.adapt_structure`; a new payoff holding arrays needs one too.
- Rates: `ZeroRate(r)` scalar or `ZeroRateCurve`/`ImpliedZeroRate` built on `Curve` (alias `CurveType`, [curves.jl](src/utils/curves.jl)): a type of the package with two sorted vectors, `times` and `values`, linear between the nodes and flat outside (`value_at`). Curves combine with `+`, `-`, `*`, `/` with numbers and with each other on the union of their times (no external dictionary type: arithmetic on one was type piracy, #42). Code uses `integral(r, T)` / `incremental_integral` so the same path works for scalar and curve rates.

### Payoffs ([src/options/](src/options/))

Hierarchy in [payoff_base.jl](src/options/payoff_base.jl): `SingleNamePayoff` → `EuropeanPayoff`, `AmericanPayoff`, `BermudanPayoff`, `PathDependentPayoff` (→ `BarrierPayoff`, `AsianPayoff`); separately `BasketPayoff`. Generic engines in [src/options/engines/](src/options/engines/) implement `payoff(S, ...)` per family (American/Bermudan use Longstaff–Schwartz); a concrete payoff usually only defines its struct and a scalar `payout(ST, payoff)` (European) or path-based equivalent. Options made of positive numeric fields plus `isCall` are defined with `@option_type Name Family T K ...` in [payoff_base.jl](src/options/payoff_base.jl), and `vanilla_payout` holds the call/put formula; barrier options only define `is_alive(S, payoff)`. American and Bermudan share one Longstaff-Schwartz `payoff` and differ in `exercise_columns`; its 3×3 normal equations are solved on scalars by `solve_3x3` (Gaussian elimination with partial pivoting, as LAPACK): no allocations, and faster than `lu!` in every backend, Enzyme included, which has no rule for LAPACK's `getrf`. The regression is on the return from the spot, `(S - S0) / S0` (`lsm_regressor`), not on the price: on prices around 100 the columns 1, S, S² are almost collinear and the fitted continuation values lost two digits (#35, tested against a `BigFloat` fit in `test/engines.jl`); the sums and the elimination use `muladd`, the continuation value `evalpoly`. Payoffs must expose `.T` (used by `maturity`).

Smiles ([src/options/smiles_base.jl](src/options/smiles_base.jl), `EuropeanOptionSmile`, `BinaryEuropeanOptionSmile`, aliases of `OptionSmile`) are not `AbstractPayoff`s (a smile has one value per strike, so it stays out of `Position`), but have the same interface: `payout(ST, smile)` returns one value per strike, `payoff(S, smile, ...)` a matrix with one row per path and one column per strike, and `pricer`, `variance`, `confinter`, `delta`, `rho` one value per strike, without building the options. A smile and its single option share one formula, `payout_strike(ST, K, isCall, OptionType)`, so a new smile family needs only that method and its struct. Strikes and flags are `AbstractVector`s, and `smile_payoffs` is one broadcast of the values at maturity against `transpose(K)` on any device (no `reshape`); one-dimensional inputs are typed `AbstractVector` throughout (Bermudan dates, curve nodes, `set_parameters!`). In `CudaMode` all the strikes run in one kernel and these metrics return a `CuArray` (`confinter` a `Vector` of tuples); `pricer` and `payoff_moments` move the smile to the mode first (`to_mode`), so a host smile is priced on the GPU and a GPU smile on the host without errors. `smile_to_options` remains for what really needs the options (VibratoMonteCarlo's `singular_points`).

### Contract algebra / portfolios ([src/operations/](src/operations/))

- `Position{P, W}` holds its legs as two vectors, `payoffs::Vector{P}` and `weights::Vector{W}`, in insertion order (a position has few legs, so a payoff is found by `isequal`; an `OrderedDict` there was untypable for Enzyme, #40): payoffs combine with `+`, `-`, `*`, `/` through `combine` (payoff type `typejoin`, weight `promote_type`, `foreach` over the tuple of terms); a weight that is exactly zero (`is_null_weight`: integers, floats, rationals only, never dual numbers) drops its leg. Payoffs compare by structure: `hash`/`isequal`/`==` on `Base.typename` and `parameter_values`, so `K = 100` and `K = 100.0` are one leg.
- A position is priced path by path: `payoff(S, position, ...)` is the weighted sum of its legs' discounted payoffs, so `pricer`, `variance`, `confinter`, `delta`, `rho` take positions; `maturity(position)` starts its fold from a zero of the maturity type (an `Int8` in the fold made a union that Enzyme rejects).
- `MarketDataSet` and `Portfolio` wrap `OrderedDict{UnderlyingKey, ...}`, keys a `String` or a tuple of names for multivariate models and baskets (`("ENI", "AAPL") → copula`); `"A_B"` is deprecated (`underlying_key`). `pricer(::MarketDataSet, ..., ::Portfolio)` seeds once and simulates the models in order from one stream (`@set mcConfig.reseed = false`), so the underlyings are independent (#38); `positions_on` gives each multivariate model its positions per subset of components (`compute_indices(n)`: the nonempty subsets of `1:n` by size, then lexicographically). Printing is `show(io, ::MIME"text/plain", ...)`, no `display` override.

### Automatic differentiation conventions

AD (ForwardDiff, ReverseDiff, Mooncake, Enzyme, DualNumbers, TaylorSeries, HyperDualNumbers) is a first-class concern, which drives several pervasive patterns — follow them in new code:
- Structs are fully parametric on each numeric field and carry a promoted `numtype` parameter (e.g. `KouProcess{num, num1, …, numtype}`), computed in the inner constructor via `promote_type` or `zero(num1) + zero(num2)`.
- Output element types come from `price_type(inputs...)` ([src/output_type.jl](src/output_type.jl)), a number type (`typeof(zero(T1) + zero(T2) + …)` over the inputs' `number_type`s), used to allocate path/price arrays of the right (possibly Dual) type — never hardcode `Float64` buffers. It works on types only (no array traversal, so GPU arrays and ReverseDiff tapes are untouched); `BaseProcess{T}` and `AbstractPayoff{T}` carry their type, anything else needs a `number_type` method (there is no catch-all, so a missing one is a `MethodError`, not a silently ignored input). `get_matrix_type`/`get_array_type(mcConfig, process, T)` take that type.
- Assertions/validation and non-differentiable bookkeeping are wrapped in `ChainRulesCore.@ignore_derivatives` (needed for Zygote/ReverseDiff compatibility).
- No Float64 literals or types in package code (`0.5 * x`, `max.(x, 0.0)`, `d = 0.0` promote Float32 inputs to Float64): use integers, `zero(x)`, `one(x)`, `x / 2`. `test/metrics/statistics.jl` checks that Float32 inputs give Float32 prices.
- Floors must stay away from the subnormal numbers (`eps(0.0)` is one, and arithmetic on subnormals is ~100× slower): Heston floors its variance at `variance_floor(zero) = zero + sqrt(floatmin(eps(zero)))`.
- Mooncake works with any RNG: `FinancialMonteCarloMooncakeExt` makes `set_seed!` a zero-derivative primitive, because seeding the default `MersenneTwister` ends in a foreigncall into libdSFMT. Mooncake's first rule build is slow (~2 min in `test/ad/mooncake.jl`).
- Enzyme (forward and reverse): `FinancialMonteCarloEnzymeCoreExt` (weak dependency EnzymeCore, loaded by Enzyme) marks `set_seed!` and the scalar draws `draw_uniform`/`draw_exponential`/`draw_normal` ([paths.jl](src/models/paths.jl)) inactive, since Enzyme fails on the scalar draws of the `MersenneTwister`: draw scalars through them, never `rand(host_rng(...))`. Loop over `path_sides` with `foreach` (a `for` over the tuple gives an `EnzymeNoTypeError`), and keep LAPACK out of the differentiated code. Variance Gamma does not work (the Gamma sampler of Distributions draws scalar exponentials). On Julia 1.10 a position built inside the differentiated function needs `Enzyme.set_runtime_activity(Reverse)` (its constant integer weights sit next to active payoffs); Julia 1.11 and later do not. Allocation tests measure inside a function (`f(x) = @allocated g(x)`): at the top level Julia 1.10 counts the boxing of the result.

### Optional dependencies

[ext/](ext/) holds the package extensions, one per weak dependency: CUDA (GPU models, in [ext/FinancialMonteCarloCUDAExt/](ext/FinancialMonteCarloCUDAExt/)), DifferentialEquations (price directly from an SDE), DualNumbers, EnzymeCore, Mooncake, TaylorSeries, VectorizedRNG, Distributed. Code there only exists when the user has loaded the corresponding package. An extension can only add methods, so types and functions it needs must be defined in `src/` (`CudaMode`, `MultiProcess`, `MonteCarloDiffEqModel`, `HostAdaptor`, `default_nbatches`, `run_on_batches`, `array_type`, `mode_adaptor`, `seed_device!`), and it imports the unexported names of the package explicitly. The CUDA extension changes behaviour only through the mode (`array_type`, `mode_adaptor`, `seed_device!`) and the array type (`CuArray`/`CuMatrix` methods), never through the configuration type. CUDA is not a test dependency, so nothing runs that extension on CI.

## Interface shared with the other packages

FinancialFFT, VibratoMonteCarlo and FinancialCalibrator (in `~/.julia/dev`) extend `FinancialMonteCarlo.pricer` with their own `AbstractMethod`s. The contract: `pricer(model, zeroCurve, method, payoff)` with four arguments, the mode inside the method (read with `execution_mode`, changed with `with_mode`), device placement through `to_mode`/`mode_adaptor` and dispatch on `AbstractCudaMode`. A change of these names is a breaking change for those packages too.

## Extending (see docs/src/extending.md for runnable examples)

- New model: struct subtyping the appropriate abstract type (engine vs. process matters), holding its parameters (+ `underlying` for processes), and a `simulate!` method; copy the signature from an existing model.
- New payoff: struct subtyping the right payoff family with its parameters (`T`, `K`, barriers…) and a `payoff`/`payout` method; all metrics then work automatically.
- `get_parameters` / `set_parameters!` exist on models and payoffs and are tested, inferred included ([test/models/parameters.jl](test/models/parameters.jl)) — keep them consistent when adding fields. They read and rebuild the fields through `parameter_names`/`parameter_values` ([src/io_utils.jl](src/io_utils.jl)), generated so the names are constants: never loop over `fieldnames` with `PropertyLens{symb}` or `getproperty(x, symb)`, a runtime Symbol makes the result `Any`. A model whose parameters are not one field each (LogNormalMixture) needs its own `set_parameters!`.
- Results must be inferred also for an array of payoffs of different types (abstract element type): take the element type from `price_type` and annotate the comprehension (`pricer`, `payoff_moments`). A `Position`, a `Portfolio` or a `Vector{AbstractPayoff}` holds abstract types, so their prices are not inferred.

## Documentation

The design is shared by the docs of all the rcalxrc08 packages (identical copies in every repository, edit them together): `docs/src/assets/style.css` (Dark Quant palette), `dark.js` (dark theme by default), `figures.css`, `site.css` (page headers with kicker, the first paragraph as lede, chips and a strip of section links; numbered sections; figure, maths, table, docstring and footer cards) and `home.css` (landing page: hero, stat tiles, feature cards, flow diagram, results behind radio-button tabs, quick start steps). `index.md` is that landing page; `docs/figures.jl` (run by `make.jl` in a module of its own) simulates its figures at build time, a light and a dark SVG each in `docs/src/assets/figures/` (git-ignored): paths, convergence with `confinter`, terminal densities, and the hero drawn by `docs/surface.jl`, the density of the Heston log-return over time. The tutorials are generated by Literate and have no page header.

Keep the docs in step with the code in the same commit: the manual pages in `docs/src/manual/` (configuration, models, payoffs, autodiff, performance), `docs/src/api.md` (every exported or documented name), `docs/src/extending.md`, the Literate tutorials in `docs/literate/`, the docstrings and this file. Measured numbers in the docs (performance and autodiff tables) say on which machine they were taken; refresh them when a change moves them.

## Releasing

1. Bump `version` in `Project.toml` on master (breaking changes: minor bump while 0.x), wait for CI to be green on Julia 1.10 and 1.x.
2. Register by commenting on the release commit on GitHub, with release notes (sections: Breaking changes, New features, Performance, Bug fixes):
   `gh api repos/rcalxrc08/FinancialMonteCarlo.jl/commits/<sha>/comments -F body=@notes.md`, where `notes.md` starts with `@JuliaRegistrator register` then `Release notes:`. JuliaRegistrator answers on the commit with the General PR; AutoMerge merges it in ~20 minutes.
3. **Tag by hand**, TagBot usually does not: `git tag -a vX.Y.Z -m "FinancialMonteCarlo vX.Y.Z" <sha>`, `git push GitHub vX.Y.Z`, then `gh release create vX.Y.Z -R rcalxrc08/FinancialMonteCarlo.jl --title vX.Y.Z --notes-file <notes without the registrator lines> --verify-tag`.
4. Check that `https://rcalxrc08.github.io/FinancialMonteCarlo.jl/stable/` serves the new version (the tag run starts a run on master, which builds `stable/` from the highest `v*` tag; see below).

Push to the `GitHub` remote; `origin` (GitLab) has not been pushed since 0.5.2 and is many commits behind.

### Docs deployment

- GitHub Pages is set to **build from GitHub Actions** (`build_type: workflow`, `gh api repos/rcalxrc08/FinancialMonteCarlo.jl/pages`); the `source.branch: gh-pages` in that answer is a leftover and is not served. There is no `gh-pages` branch any more (#52): the `docs` job of [Documentation.yml](.github/workflows/Documentation.yml) builds the manual of master into `docs/build` (the site root), checks out the highest `v*` tag (`git ls-remote --tags | sort -V`) into `release/` and builds its manual into `docs/build/stable` (with `GITHUB_REPOSITORY` unset, so the `deploydocs` of the `make.jl` of releases up to 0.7 does nothing), adds `docs/redirect/404.html`, and uploads `docs/build` as the Pages artifact, which the `pages` job deploys with `actions/deploy-pages`. `docs/make.jl` has no `deploydocs`, and its `canonical` is the site root.
- Old URLs of the versioned deploydocs site: `dev/...` and `vX.Y.Z/...` no longer exist, and `404.html` redirects them with JavaScript to the same page at the root (`dev/`) or in `stable/` (a version folder). `stable/...` keeps working. Docs of releases older than the latest are not published: build them locally from their tag.
- A tag push runs only the `publish-tag` job, which starts the workflow on master (Pages refuses deployments from a tag; a `workflow_dispatch` started with `GITHUB_TOKEN` does start a run). Pull requests build the manual of master only, without deploying.
- The CI and Documentation runs of 54ecd0b (2026-11-11) stayed `queued` for hours with Actions reported operational, and `gh run cancel` answered that they were completed: the next push started fresh runs normally.

## Caveats (from README)

- First- and second-order AD "works" for every payoff, including ones where it is meaningless (binary, barrier); second-order results are mostly useless.
- No multicurrency support.
