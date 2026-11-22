# CLAUDE.md

**FinancialMonteCarlo.jl** prices equity derivatives and their risk by Monte Carlo, single-name and multivariate
(Gaussian copula), through `pricer(model, zeroCurve, method, payoff)`, the interface FinancialFFT, VibratoMonteCarlo and
FinancialCalibrator extend with their own methods. Every model and payoff is differentiable by every Julia AD back-end
and runs serially, on threads, on processes and on CUDA. Closed forms in tests and docs come from FinancialToolbox.
Registered in General. Goals and non-goals: `references/mission.md`; why things are the way they are:
`references/decisions.md`.

This workspace follows the Interpretable Context Methodology (layered context, arXiv:2603.16021): load only what the
task needs.

| Layer | Where | What |
|:--|:--|:--|
| 0 | this file | identity, map, rules that hold everywhere |
| 1 | `CONTEXT.md` | which stage or folder a task belongs to, shared references |
| 2 | `src/0N_*/CONTEXT.md`, `ext/`, `test/`, `docs/`, `benchmark/`, `etc/` | a contract per stage: inputs, process, outputs, runtime checks |
| 3 | `references/*.md` | stable facts, decisions and measurements, loaded when a contract lists them |
| 4 | `state/STATUS.md` | what is in flight: current focus, open issues, known problems |

Startup: read `state/STATUS.md`, then route the task through `CONTEXT.md`. Runtime truth comes from a Kaimon session
for this project (`kaimon.toml` pins the Julia); each contract's Verify table lists the checks to run after a change.

Map: `src/FinancialMonteCarlo.jl` includes the stages in order, `01_utils` → `02_models` → `03_options` →
`04_operations` → `05_metrics` → `06_batches` (each keeps its original aggregator file and folder, relative includes
unchanged); `ext/` one extension per weak dependency; `test/` the suite (`test/cuda/` local GPU tests); `docs/` the
Documenter + Literate site; `benchmark/` the suite; `etc/` the logo.

Rules for every change (details in `references/conventions.md` and `references/autodiff.md`):

- Every method carries its execution mode; the four-argument `pricer` is shared with three packages (breaking if
  changed).
- One `simulate!` body per model on the hooks of `paths.jl`; inner loops as kernels with closed-form adjoints and
  Mooncake/Enzyme rules.
- Parametric structs with a promoted `numtype`; output types from `price_type`; no Float64 literals or buffers.
- Validation by `check_domain`/`check_dimension`/`check_argument`, never `@assert`, inside one named function per
  constructor marked `@non_differentiable` (no `@ignore_derivatives` block, D-013); `@muladd` for every `a * b + c`.
- Reuse the author's packages (closed forms from FinancialToolbox) and established libraries.
- Format with `.JuliaFormatter.toml`; update the CONTEXT.md, references, docs and `state/STATUS.md` a change touches in
  the same commit (a new rule also gets an entry in `references/decisions.md`); push to the `GitHub` remote.
