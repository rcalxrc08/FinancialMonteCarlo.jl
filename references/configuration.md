# Stage 01: configuration, curves, underlyings

_Methods, modes, samplings, the configuration, data placement, curves._ (Layer 3 reference.)

## Three concepts (`mc_mode.jl`, `mc_config.jl`)

- **Method** (`AbstractMethod`): what `pricer(model, underlying, zeroCurve, method, payoff)` dispatches on.
  `MonteCarloConfiguration <: AbstractMonteCarloConfiguration <: AbstractMethod`; FinancialFFT, VibratoMonteCarlo and
  FinancialCalibrator subtype it too. **Every method carries its execution mode**: `execution_mode(method)` reads it,
  `with_mode(method, mode)` returns a copy in another mode (generic `@set method.mode`; `MonteCarloConfiguration`
  rebuilds its random state). Never add a trailing `mode` argument to a metric.
- **Mode** (`AbstractMode`): where and how a method runs, with no seed and no generator, so methods without randomness
  use the same objects. `SerialMode()`, `MultiThreading(nbatches)`, `MultiProcess(nbatches)` (batched modes,
  `AbstractBatchedMode`), `CudaMode()` (`AbstractCudaMode`, defined in `src/` so other packages can dispatch on it
  without CUDA).
- **Sampling** (`AbstractSampling`): `StandardMC`, `AntitheticMC`, `SobolMC`, or `ControlVariates` (defined in
  `payoff_base.jl`, handled in `pricer_cv.jl`, `metrics.md`). The aliases `SobolConfig` and `BatchedConfig` are defined
  in `mc_config.jl`.

## `MonteCarloConfiguration`

- `MonteCarloConfiguration(Nsim, Nstep; sampling = StandardMC(), mode = SerialMode(), seed = 0, rng = MersenneTwister(),
  reseed = true)` is the only public constructor (fields `Nsim, Nstep, sampling, mode, rng, reseed`; the default
  6-argument constructor is internal, used by `batch_configs` and `@set`).
- It runs every compatibility check once, in `check_configuration` (`check_sampling`, `check_sampling_mode`,
  `check_reseed`, `check_mode`), so an unsupported combination fails at construction; add new restrictions there as a
  dispatch method, not in the metrics.
- Its `rng` field is the random state built by `random_state(mode, seed, rng)`: one `SeededRNG(seed, rng)`, or for a
  batched mode a `Vector{SeededRNG}` (batch `i` gets a copy of `rng` with seed `seed + i`, or `seed[i]` when `seed` is a
  vector). `set_seed!` reseeds them (if `reseed`) and calls `seed_device!(mode, seed)` (CURAND in the CUDA extension).
- The 0.5 names (`SobolMode`, `BaseMode`, `AbstractMonteCarloMethod`, positional constructors) were removed in 0.9 (#77).

## Data placement

The mode decides where the computation runs and where results live. Inputs holding arrays are moved to the mode once
when a metric starts, with `to_mode(x, mode)` (`Adapt.adapt(mode_adaptor(mode), x)`; `HostAdaptor` for the CPU modes,
`CuArray` for CUDA). Smiles define `Adapt.adapt_structure`; a new payoff holding arrays needs one too.

## Rates and curves (`zero_rates.jl`, `curves.jl`, `underlying.jl`)

- One type per concept since 0.10 (D-014): `ZeroRate{R}(r)` and `Underlying{S, D}(S0, d)`, the rate and the dividend a
  number or an `AbstractCurve`, a type parameter; methods that differ dispatch on it (`ZeroRate{<:Number}`,
  `ZeroRate{<:Curve}`, e.g. `bump_rate`). A curve is a `Curve{N, T1, T2}` (or an `ImpliedCurve` from quoted zero
  rates): two sorted `SVector`s (StaticArrays, a hard dependency since 0.9), `times` and `values`, linear between the
  nodes and flat outside (`value_at`). The nodes are static so that a curve and an underlying or rate holding one are
  isbits (#79): Enzyme needs no runtime activity, GPU kernels capture them.
- The algorithms are the generic `AbstractVector` ones (`searchsortedlast`, early returns): a tuple prototype with
  branch-free loops was slower at large N and compiled for tens of seconds (measurements on #79). N is in the type:
  built from vectors (`Curve(r, T)`, `ImpliedCurve`) the type is known only at run time (one dispatch where the curve is
  built, e.g. FinancialCalibrator's `adapter`), from tuples it is inferred.
- Curves combine with `+`, `-`, `*`, `/` with numbers (`map` on the values); between two curves only `+` and `-`, which
  return a lazy `CurveDifference` (`integral` is the difference of the integrals, exact; the drift `r - d` when both
  are curves), no `*` or `/` (#81). `AbstractCurve` is the supertype of both, and what the engines
  (`BrownianMotionVec`, …) and `step_drifts` take; `Underlying` and `ZeroRate` hold a number or a `Curve`. No external
  dictionary type: arithmetic on one was type piracy (#42).
- Code uses `integral(r, T)` / `incremental_integral` so the same path works for scalar and curve rates.
  `integral(::Curve, ::Number)` has a ChainRules `rrule` (`integral_partials`, `autodiff.md`).

## Validation (`checks.jl`)

Input validation never uses `@assert` (#71): `check_domain(cond, value, msg)` (`DomainError`), `check_dimension(cond,
msg)` (`DimensionMismatch`) and `check_argument(cond, msg)` (`ArgumentError`), so callers (FinancialCalibrator's loss)
can catch invalid parameters without swallowing bugs.
