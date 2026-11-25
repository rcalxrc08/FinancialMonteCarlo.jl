# Stage 02: models and simulation

_The model hierarchy, the shared `simulate!` body, Sobol, the device draws._ (Layer 3 reference; the per-model
catalogue is `docs/src/manual/bestiary.md`.)

## Hierarchy (`base_models.jl`)

- `AbstractMonteCarloEngine`: raw engines with no notion of rates/dividends (`BrownianMotion`,
  `GeometricBrownianMotion`, `SubordinatedBrownianMotion`); `simulate!(X, engine, mcConfig, T)`.
- `AbstractMonteCarloProcess`: financial models holding an `underlying::AbstractUnderlying` (S0, dividend);
  `simulate!(X, process, zeroCurve, mcConfig, T)`. Subtypes: `ItoProcess` (BlackScholes, Heston, LogNormalMixture…),
  `FiniteActivityProcess` (Kou, Merton), `InfiniteActivityProcess` (VarianceGamma, NIG), and
  `VectorialMonteCarloProcess`/`NDimensionalMonteCarloProcess` (copula-based `GaussianCopulaNVariateProcess` in
  `nvariate*.jl`).
- Wrappers: `BachelierProcess` is a `BrownianMotion` (zero drift) times the carry factors `carry_factors(μ, mcConfig, T)`
  (`e^{∫₀ᵗ (r - d)}` on the grid, a row on the device of the mode, `paths.jl`); `ShiftedProcess(model, α)` wraps any
  `ScalarMonteCarloProcess`: it simulates `model` from `Underlying((1 - α) S0, d)`, built from the underlying it is
  given, then `add_shift!`. `ShiftedLogNormalMixture` was removed: it is `ShiftedProcess(LogNormalMixture(η, λ), α)`.
- `LogNormalMixture` is the Brigo–Mercurio local volatility model (#56): its `η` and `λ` are tuples (`NTuple{N}`,
  `NTuple{N-1}`), the number of components is in the type, `lnm_local_variance` works on one path with the tuples, so the
  same scalar function runs inside one broadcast on CPU and GPU (the tuples passed as `(η,)`) and in the Sobol loop;
  Euler on the log price, the local variance from log-sum-exp weights, its limit at `t = 0`.
- Models are thin: they compute a risk-neutral drift and delegate to an engine, then exponentiate. `fa_levy.jl`
  implements `simulate!` for all finite-activity Lévy models (Kou, Merton provide `compute_jump_size` and
  `compute_drift`); `ia_levy.jl` for infinite-activity ones (VarianceGamma, NIG provide `compute_drift` and
  `compute_subordinator_law`, each in two forms: one reading the model's fields, delegating to one taking the parameters
  explicitly, e.g. `compute_drift(model, σ, θ, κ)`, `compute_subordinator_law(model, κ, dt)`, which can be called with
  other values and numeric types, e.g. to build a subordinator without Duals).

## One `simulate!` body per model (`paths.jl`)

A model has one `simulate!` body for antithetic variates, constant and curve rates, and CPU and GPU, built on the hooks
of `paths.jl`; do not add a copy per variant:

- `step_drifts(μ, dt, Nstep)` / `drift_of_step(drifts, j)`: a constant rate stays one `Number` and only a curve is a
  vector, so flat simulations allocate nothing for the drift.
- `path_sides(X, mcConfig)`: all the paths or, with `AntitheticMC`, the two halves with the sign of their normal draws.
- `draw_buffer`, `draw_normals!`, `draw_subordinator!`, `get_matrix_type`: buffers and draws on the device of the mode.
  The arrays come from `array_type(mode)`; a device adds its method of `array_type` and, for the arrays of its own type
  (`Z::CuArray`), of `draw_normals!` and `draw_subordinator!(dt_s, law, mcConfig)`. Models draw on the host from
  `host_rng(mcConfig)`.
- What CUDA keeps off the host: the `Gamma` and `InverseGaussian` subordinators are drawn by a kernel (`device_draw`,
  Marsaglia-Tsang and Michael-Schucany-Haas; any other law is sampled on the host by Distributions and copied), whose
  generator is seeded at every launch from `host_rng`, because `CUDA.seed!` and `Random.seed!` do not reach the generator
  of the kernels. The jumps of Kou and Merton are drawn on the host into a `JumpList` (the samplings write through
  `add_jump!(sink, i, jump_idx, size)`, the sink being the matrix or the list) and added by one kernel with a thread per
  element (contiguous writes; a thread per path writing along its row was 4 times slower).
- Floors stay away from the subnormal numbers (`eps(0.0)` is one, and arithmetic on subnormals is ~100× slower): Heston
  floors its variance at `variance_floor(zero) = zero + sqrt(floatmin(eps(zero)))`.

## Sobol

Sobol is the one separate body (`*_sobol.jl`, dispatching on `SobolConfig = MonteCarloConfiguration{<:SobolMC}`),
because it draws path by path (the first `2^m - 1` points are skipped by `sobol_skip!`, which sets the state of the
sequence in closed form: never call `Sobol.skip`, which draws every skipped point); it fills each path inside
`on_sobol_paths`, a function barrier, because the type of a `SobolSeq` holds its dimension, known only at run time (the
one intended runtime dispatch JET reports).

## Jumps from uniforms

Kou and Merton (`FiniteActivityProcess`) draw their jumps from uniforms by inversion with `AntitheticMC` (its table of
counts is in `Float64`, the type of `rand(rng)`, without a draw for the type, #54) and `SobolMC` (`fa_levy.jl`): the
count from the Poisson table `JumpUniforms` (`jump_count`), then a time and a size (`jump_quantile`, one method per
model) per jump; antithetic paths take `1 - u`, and a Sobol point has `Nstep + 1 + 2M` coordinates
(`on_sobol_paths(...; normals = Nstep)` leaves the jump coordinates uniform), beyond `M` jumps the rng. VibratoMonteCarlo
reuses these helpers. Under AD the jumps are recorded first (`autodiff.md`, "Kernels with rules").

## Parameters (`src/04_operations/io_utils.jl`)

`get_parameters` / `set_parameters!` exist on models and payoffs and are tested, inferred included
(`test/models/parameters.jl`); keep them consistent when adding fields. They read and rebuild the fields through
`parameter_names`/`parameter_values`, generated so the names are constants: never loop over `fieldnames` with
`PropertyLens{symb}` or `getproperty(x, symb)`, a runtime Symbol makes the result `Any`. A model whose parameters are
not one field each (LogNormalMixture: two tuples, rebuilt with `ntuple(…, Val(N))` so it is inferred) needs its own
`set_parameters!`. `parameter_count(T)` is the length of `get_parameters` from the type alone (a constant: it sizes the
ForwardDiff chunk of FinancialCalibrator, FinancialCalibrator#36); an array parameter throws `ArgumentError`. `parameter_fields` excludes only the tuple of models of a copula (`Tuple{Vararg{BaseProcess}}`), so
tuples of numbers are parameters (`flat_parameter(::Tuple)`). `get_parameters(::ShiftedProcess)` has its own method
(the parameters of the model, then `α`): through the generic one it would call itself on the inner model and inference
gives up.

## Adding a model

A struct subtyping the appropriate abstract type (engine vs process matters), holding its parameters (+ `underlying`
for processes), and a `simulate!` method (copy the signature from an existing model); its step as a kernel with an
adjoint and two rules (`autodiff.md`); a card in `docs/src/manual/bestiary.md`; models easy to add are tracked in #62.
