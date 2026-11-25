# Decisions

_Architecture decisions, oldest first. Append; a reversed decision gets a new entry that says "supersedes D-NNN"._
(Layer 3 reference.) Read before a refactor: a structure that looks accidental may be an entry here.

## D-001 (0.6, FinancialFFT#55) The execution mode is a field of the method

- **Decision:** methods carry their mode (`execution_mode`, `with_mode`); no trailing `mode` argument.
- **Why:** one dispatch path for every package extending `pricer`.
- **Consequences:** the shared interface of `public_api.md`.
- **Status:** active

## D-002 One `simulate!` body per model

- **Decision:** antithetic, curves, CPU and GPU through the hooks of `paths.jl`; Sobol is the only separate body.
- **Why:** copies per variant drifted; the hooks keep flat simulations allocation-free.
- **Consequences:** `models.md`.
- **Status:** active

## D-003 (#71) Explicit exceptions, never `@assert`

- **Decision:** `check_domain`, `check_dimension`, `check_argument`.
- **Why:** FinancialCalibrator's loss catches invalid parameters without swallowing bugs.
- **Status:** active

## D-004 Kernels with closed-form adjoints and Mooncake/Enzyme rules

- **Decision:** inner loops as kernels (`bm_step!`, `heston_step!`, …) with adjoints in `adjoints.jl` and rules in the
  extensions; jumps recorded before their sizes.
- **Why:** without them Mooncake reverse on Heston and the mixture was 3-7 times the price.
- **Consequences:** a new model needs a kernel, an adjoint, two rules and AD tests (`autodiff.md`).
- **Status:** active

## D-005 (#35) Longstaff–Schwartz on the return, `solve_3x3`

- **Decision:** regress on `(S - S0) / S0`; solve the 3×3 system on scalars.
- **Why:** two digits lost on prices around 100; LAPACK has no Enzyme rule and allocates.
- **Status:** active

## D-006 (#40) `Position` as two vectors

- **Decision:** legs as `payoffs::Vector{P}`, `weights::Vector{W}`, not an `OrderedDict`.
- **Why:** the dictionary was untypable for Enzyme.
- **Status:** active

## D-007 (#76) Variance of independent samples

- **Decision:** antithetic pairs are one sample; `variance` per path.
- **Why:** treating correlated pairs as independent gave intervals 30% too narrow.
- **Status:** active

## D-008 (#79, #81; 0.9.0) `SVector` curves, no curve products

- **Decision:** curve nodes are `SVector`s (StaticArrays hard dependency); between curves only `+`/`-` (lazy
  `CurveDifference`).
- **Why:** isbits curves and models (Enzyme without runtime activity, GPU capture); the union of times had a length
  known only from the data.
- **Consequences:** FinancialFFT and FinancialCalibrator call `integral` per curve.
- **Status:** active

## D-009 (#77; 0.9.0) Removal of the 0.5 names

- **Decision:** `SobolMode`, `BaseMode`, `AbstractMonteCarloMethod`, positional constructors removed;
  `MonteCarloConfiguration` keyword constructor only.
- **Status:** active

## D-010 (2026-11-19, user) Zygote scope: the structs, not the simulations

- **Decision:** constructors, `set_parameters!`, curves, `integral`, underlyings are Zygote-differentiable (two rrules);
  simulations are not targeted.
- **Why:** FinancialCalibrator differentiates `vol_errors` through them (#80).
- **Status:** active

## D-011 (2026-11-22) Layered context (ICM)

- **Decision:** the repository context follows the Interpretable Context Methodology: `CLAUDE.md` (L0), `CONTEXT.md`
  routing (L1), stage contracts (L2), `references/` (L3), `state/` (L4); `test/context.jl` keeps it in sync (#84).
- **Why:** the 5348-word `CLAUDE.md` loaded everything for every task; the layout worked in FinancialCalibrator and
  VibratoMonteCarlo.
- **Consequences:** the top-level source files and their folders moved into stage folders (relative includes unchanged).
- **Status:** active

## D-012 (2026-11-22, #64) The SDE extension is triggered by StochasticDiffEq

- **Decision:** `MonteCarloDiffEqModel` is simulated by `FinancialMonteCarloStochasticDiffEqExt`, a weak dependency on
  StochasticDiffEq (6, 7), not on DifferentialEquations.
- **Why:** DifferentialEquations 8 contains only SciMLBase and OrdinaryDiffEq: no SDE solver (`SOSRI`), no jumps; widening
  its compat would load an extension that cannot work. StochasticDiffEq has the solver and reexports the problem types.
- **Consequences:** users of DifferentialEquations 7 get the extension as before (it loads StochasticDiffEq); with
  DifferentialEquations 8 they load StochasticDiffEq (and JumpProcesses for jumps) themselves. Release note needed.
- **Status:** active

## D-013 (2026-11-23, #85) Validations are `@non_differentiable` functions, not `@ignore_derivatives` blocks

- **Decision:** every validation and every non-differentiable bookkeeping (the `numtype` zero, the call/put sign, the
  legs of a position) is a named function marked `ChainRulesCore.@non_differentiable`, called directly with the raw
  inputs, one per constructor; the throws are `@noinline` (`throw_domain_error` & co., as FinancialFFT, Enzyme#4053).
- **Why:** an `@ignore_derivatives` block is a closure that Zygote traces: ~2.5 μs each with Zygote 0.7.19 (~11 μs with
  0.7.18), per element in a broadcast payout; the function costs nothing measurable. FinancialCalibrator rebuilds the
  model at every iteration. Enzyme, Mooncake and ForwardDiff are indifferent.
- **Consequences:** the exception types are unchanged; `check_domain`, `check_dimension`, `check_argument` are
  non-differentiable too, for the callers of a single check.
- **Status:** active

## D-014 (2026-11-25, #72; 0.10.0) The underlying is an argument, not a field of the models; one type per market object

- **Decision:** models hold only their parameters; every entry point takes the `Underlying` right after the model,
  `pricer(model, underlying, zeroCurve, method, payoff)` (`simulate`, `simulate!`, `variance`, `confinter`, `delta`,
  `rho`, `distribution` alike), a tuple of underlyings for a multivariate model. A `MarketDataSet` holds
  `name → underlying` and `key → model` (`name → (model, underlying)` builds both). `Underlying{S, D}` and
  `ZeroRate{R}` are single types whose dividend and rate are a number or a curve (no `UnderlyingScalar`/`UnderlyingVec`,
  `ZeroRateCurve`, `ImpliedZeroRate`, nor vector-and-maturity constructors: `ZeroRate(Curve(r, T))`). Hard break, no
  shim: the dependent packages move in lockstep.
- **Why:** the reasons the underlying went into the models in 2019 (one spot per copula component, names as keys) are
  gone (tuples, `MarketDataSet`). As an argument it gives AD per-argument activity (delta and vega no longer share an
  object, a dual spot does not promote the model's `numtype`), one `set_parameters!` per calibration with a scalar
  dividend per maturity, bumps without `@set model.underlying`, and drops the `underlying_tangent` rules (#80, FC#44).
  One type per market object: dispatch on the field type instead of constructors returning different types.
- **Consequences:** breaking for FinancialFFT, VibratoMonteCarlo, FinancialCalibrator (0.10.0). `delta` bumps the
  underlying; in a market data set it bumps the name's underlying, shared by every model using it.
- **Status:** active
