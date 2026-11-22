# Stage 02: models and simulation

Engines and financial models, the shared simulation hooks, Sobol bodies, jumps, and the AD kernels with their adjoints.

## Inputs

| Load | Why |
|:--|:--|
| `models.jl` | includes `models/models.jl` |
| `models/models.jl` | the aggregator: the order of the files below |
| `models/base_models.jl` | the hierarchy: `BaseProcess`, `AbstractMonteCarloEngine`, `AbstractMonteCarloProcess`, `ItoProcess`, `LevyProcess`, `FiniteActivityProcess`, `InfiniteActivityProcess`, `ScalarMonteCarloProcess`, `VectorialMonteCarloProcess`, `NDimensionalMonteCarloProcess`, `simulate` |
| `models/utils.jl` | `gausscopulagen` |
| `models/paths.jl` | the hooks (`step_drifts`, `path_sides`, `draw_buffer`, `draw_normals!`, `draw_subordinator!`, scalar draws), the kernels `bm_step!`, `exp_paths!`, `affine_carry!`, `add_carry!`, `carry_factors` |
| `models/brownian_motion.jl`, `models/brownian_motion_aug.jl`, `models/brownian_motion_sobol.jl` | `BrownianMotion` (plain, vector/curve drift, Sobol) |
| `models/geometric_brownian_motion.jl`, `models/geometric_brownian_motion_aug.jl` | `GeometricBrownianMotion`, `GeometricBrownianMotionVec` |
| `models/black_scholes.jl`, `models/bachelier.jl` | Black–Scholes; Bachelier (Brownian motion times the carry factors) |
| `models/heston.jl`, `models/heston_sobol.jl` | Heston: `heston_step!`, `heston_log_step`, `heston_variance_step`, `variance_floor`; its Sobol body |
| `models/log_normal_mixture.jl` | `LogNormalMixture` (tuples `η`, `λ`), `lnm_local_variance`, `lnm_step!` |
| `models/fa_levy.jl`, `models/kou.jl`, `models/merton.jl` | finite-activity Lévy: one `simulate!`, `JumpList`, `JumpUniforms`, `record_jumps`, `add_jump_sizes!`; Kou and Merton supply `compute_jump_size`, `compute_drift`, `jump_quantile` |
| `models/subordinated_brownian_motion.jl`, `models/subordinated_brownian_motion_aug.jl`, `models/subordinated_brownian_motion_sobol.jl` | `SubordinatedBrownianMotion` (plain, curve, Sobol) |
| `models/ia_levy.jl`, `models/variance_gamma.jl`, `models/normal_inverse_gaussian.jl` | infinite-activity Lévy: one `simulate!`; VG and NIG supply `compute_drift`, `compute_subordinator_law` in two forms |
| `models/diff_eq_monte_carlo.jl` | `MonteCarloDiffEqModel` (the StochasticDiffEq extension simulates it) |
| `models/shifted_process.jl` | `ShiftedProcess(model, α)`, `add_shift!` |
| `models/adjoints.jl` | closed-form adjoints `*_adjoint!` of the kernels, used by the Mooncake and Enzyme rules |
| `models/nvariate.jl`, `models/nvariate_log.jl` | `GaussianCopulaNVariateProcess`, `GaussianCopulaNVariateLogProcess` |
| `references/models.md` | hierarchy, hooks, Sobol, jumps, parameters, adding a model |
| `references/autodiff.md` | kernels with rules, Enzyme draws, Mooncake seeding |
| `references/conventions.md` | before adding a model |

## Process

- One `simulate!` body per model (antithetic, curves, CPU, GPU) on the hooks; Sobol the only separate body, through the
  `on_sobol_paths` function barrier; never `Sobol.skip`.
- A new model: kernel, adjoint in `models/adjoints.jl`, `rrule!!`/`frule!!` (Mooncake) and
  `augmented_primal`/`reverse` (EnzymeCore), AD tests against ForwardDiff to 1e-10, a bestiary card, `get_parameters`
  consistency.
- Scalar draws through `draw_uniform`/`draw_exponential`/`draw_normal`; `foreach` over `path_sides`.
- `models/ivm/` is an untracked experiment, not included.

## Outputs

`simulate` returns the `Nsim × (Nstep+1)` path matrix on the device of the mode. Checked by `test/models/*.jl`,
`test/methods/sobol.jl`, `test/ad/*.jl`, `test/cuda/` (local).

## Verify (Kaimon)

Setup: `using FinancialMonteCarlo, Test`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.0))`,
`mc = MonteCarloConfiguration(1000, 10)`.

| Check | Expect |
|:--|:--|
| `size(simulate(bs, ZeroRate(0.01), mc, 1.0))` | `(1000, 11)` |
| `@inferred simulate(bs, ZeroRate(0.01), mc, 1.0)` | `Matrix{Float64}` |
| `eltype(simulate(BlackScholesProcess(0.2f0, Underlying(100f0, 0f0)), ZeroRate(0.01f0), mc, 1f0))` | `Float32` |
