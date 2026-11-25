# Automatic differentiation

_The conventions every new struct and kernel follows, the rules, each back-end's constraints._ (Layer 3 reference.)
AD (ForwardDiff, ReverseDiff, Mooncake, Enzyme, DualNumbers, TaylorSeries, HyperDualNumbers) is a first-class concern.

## Structs and types

- Structs are fully parametric on each numeric field and carry a promoted `numtype` parameter (e.g. `KouProcess{num,
  num1, …, numtype}`), computed in the inner constructor via `promote_type` or `zero(num1) + zero(num2)`.
- Output element types from `price_type` (`metrics.md`). No Float64 literals or types in package code (`0.5 * x`,
  `max.(x, 0.0)`, `d = 0.0` promote Float32 inputs to Float64): use integers, `zero(x)`, `one(x)`, `x / 2`.
  `test/metrics/statistics.jl` checks that Float32 inputs give Float32 prices.
- Validation and non-differentiable bookkeeping are named functions marked `ChainRulesCore.@non_differentiable`, called
  directly with the raw inputs, one per constructor (`check_heston`, `check_curve`, the `check_<Option>` the
  `@option_type` macro generates, `zero_of_types` for the `numtype`, `call_put_sign`; `checks.jl`, #85): never an
  `@ignore_derivatives` block, whose closure Zygote traces (D-013).

## Kernels with rules (the AD performance of the package)

- The inner loops of the models are functions of their own, `bm_step!`, `exp_paths!`, `affine_carry!`, `add_carry!`
  (`paths.jl`), `heston_step!` (with the scalar `heston_log_step`/`heston_variance_step`), `lnm_step!`,
  `add_jump_sizes!` (`fa_levy.jl`), and their adjoints in closed form on whole arrays are in `adjoints.jl`
  (`*_adjoint!`: accumulate into the shadows of the array inputs, `nothing` for a constant one, reset the adjoint of what
  the kernel overwrote, return the scalar derivatives).
- The Mooncake extension gives each a `rrule!!` (saving and restoring what it overwrites, as Mooncake's own in-place
  rules) and a `frule!!`; the EnzymeCore extension an `augmented_primal`/`reverse` (tape: the values the adjoint reads,
  copied when `overwritten(config)` says so; Enzyme's forward mode differentiates the kernels themselves). The rules are
  for arrays of `Float16/32/64` only; dual numbers and the GPU go through the generic broadcasts. Without them Mooncake
  records every element of every broadcast and Enzyme caches every overwritten value: Mooncake reverse on Heston and the
  lognormal mixture was 3-7 times the price.
- The jumps of Kou and Merton are drawn first as random numbers (`record_jumps`, a `JumpRecord` of grid points and
  `RawJump`/`QuantileJump` draws, inactive / zero derivative), their sizes computed after (`jump_size(model, raw)`,
  `jump_quantile`), and added by `add_jump_sizes!`: no loop of random length per path under AD.
- A new model gets the same treatment: its step as a kernel, the adjoint in `adjoints.jl`, the two rules, and a case in
  `test/ad/mooncake.jl` and `test/ad/enzyme.jl` (reverse and forward against ForwardDiff to 1e-10).
- Kernels use `@muladd` (MuladdMacro), and so do the Sobol loops (Heston's calls the scalar
  `heston_log_step`/`heston_variance_step` of its kernel), the curves and the device samplers of the CUDA extension
  (#78); a rule that recomputes a primal value (Heston's floor test) must call the same scalar function as the kernel.

## Mooncake

Works with any RNG: `FinancialMonteCarloMooncakeExt` makes `set_seed!` a zero-derivative primitive, because seeding the
default `MersenneTwister` ends in a foreigncall into libdSFMT. Mooncake's first rule build is slow (~2 min in
`test/ad/mooncake.jl`).

## Enzyme (forward and reverse)

- `FinancialMonteCarloEnzymeCoreExt` (weak dependency EnzymeCore, loaded by Enzyme) marks `set_seed!` and the scalar
  draws `draw_uniform`/`draw_exponential`/`draw_normal` (`paths.jl`) inactive, since Enzyme fails on the scalar draws of
  the `MersenneTwister`: draw scalars through them, never `rand(host_rng(...))`.
- Loop over `path_sides` with `foreach` (a `for` over the tuple gives an `EnzymeNoTypeError`), and keep LAPACK out of the
  differentiated code (`solve_3x3`).
- Variance Gamma does not work (the Gamma sampler of Distributions draws scalar exponentials).
- On Julia 1.10 a position built inside the differentiated function needs `Enzyme.set_runtime_activity(Reverse)` (its
  constant integer weights sit next to active payoffs); Julia 1.11 and later do not (#66 drops the 1.10 workarounds).
- Curve underlyings (#79, fixed in 0.9): with `Vector` nodes a model built from active parameters on an underlying on a curve
  needed `set_runtime_activity` (Enzyme 0.13.221); with `SVector` nodes it does not (tested in `test/ad/enzyme.jl`, also
  a `ShiftedProcess` of a mixture through `set_parameters!`). `EnzymeRules.inactive_type(::Type{<:Curve})` was rejected:
  it fixed only the one-level case and would zero the curve sensitivities. A closure annotated `Const` that captures
  constant heap data (FinancialCalibrator's loss cache) still makes Enzyme ask for runtime activity on any underlying:
  pass such data as `Const` arguments.
- Allocation tests measure inside a function (`f(x) = @allocated g(x)`): at the top level Julia 1.10 counts the boxing
  of the result.

## Zygote scope (user, 2026-11-19)

- The simulations are not meant for Zygote, but the structs are: model constructors, `set_parameters!`,
  `Curve`/`ImpliedCurve`, `integral`, the underlyings. FinancialCalibrator differentiates its `vol_errors` through them
  with Zygote (#80, reopened).
- Two `ChainRulesCore.rrule`s (the only ones in `src/`): `integral(::Curve, ::Number)` (closed-form partials in
  `integral_partials`, `curves.jl`: trapezoid weights, `value_at`, the times; 483 → 4 μs on the DAX curve) and
  `set_parameters!` (`io_utils.jl`: the tangent of the returned model read back in `get_parameters` order by the
  generated `parameter_tangent_values`, `LogNormalMixture` and `ShiftedProcess` methods; the underlying's tangent passed
  through; 0.56-0.67 ms → 11-19 μs). Both tested by ChainRulesTestUtils' `test_rrule` in `test/ad/chain_rules.jl` (not
  on a node, where the integral has a kink).
- #85: every `@ignore_derivatives` block costs ~11 μs under Zygote 0.7.18 (its closure, Zygote#1753), ~2.5 μs under
  0.7.19, against nothing measurable for a `@non_differentiable` function. All 111 blocks replaced (2026-11-23, Zygote
  0.7.19, i7-8550U, gradient, BenchmarkTools minimum, master against the branch in interleaved processes, two rounds):
  HestonProcess and EuropeanOption built in the function 79-134 → 25-26 μs, KouProcess 29-30 → 20.5 μs,
  `set_parameters!` on a LogNormalMixture 23-25 → 8.2-8.6 μs, the payout of a European option 62-69 → 6.8-7.0 μs,
  broadcast over 1000 values 5.7 ms → 32-37 μs (the `zero_typed` closure ran per element); ForwardDiff unchanged
  within noise (1.3-2.4 μs). Script: `benchmark/bench_validation_85.jl`.
- Zygote traps: `oftype` has a branch and an uninferred pullback (FluxML/Zygote.jl#1751, still open: use
  `promote`/`convert`). Fixed in Zygote 0.7.19 (2026-11-22): `gradient`/`pullback` compiled for `f::Function` (#1752),
  `Core._typeof_captured_variable` of every closure capture traced in reverse mode (#1753), the `Vector{Any}` stacks of
  branch-only code (#1736, #1754: a ternary on a flag 18-21 → 2.8-2.9 μs). Benchmark Zygote with interpolated non-constant arguments
  and no `$(args)...` splat.

## TaylorSeries

Needs no extension: its `isless` makes the `max`/`min` of Base work, and payoffs take the zero of the paths (#53).

## GPU

Only ForwardDiff works in `CudaMode` (#68): Enzyme forward/reverse and Mooncake reverse are `@test_broken` in
`test/cuda/test_ad_cuda.jl`, Mooncake forward crashes and is skipped. Enzyme forward one direction at a time runs on
CUDA < 6.5 but returns wrong derivatives and primal (a host-to-device `copyto!` loses its shadow, JuliaGPU/CUDA.jl#3332);
CUDA ≥ 6.4.2 needed LLVM.jl 10, which Enzyme did not allow (EnzymeAD/Enzyme.jl#3741; lifted in Enzyme 0.13.221, then
CUDA 6.5.1 asserts, EnzymeAD/Enzyme.jl#4094).
