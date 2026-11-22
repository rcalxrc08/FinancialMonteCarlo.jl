# Stage 05: metrics

The pricing pipeline and every metric built on it, and the output number types.

## Inputs

| Load | Why |
|:--|:--|
| `metrics.jl` | the aggregator |
| `metrics/seed.jl` | `set_seed!`, `reseed!`, `seed_device!` |
| `metrics/utils.jl` | `compute_indices` (included by `metrics/pricer.jl`) |
| `metrics/pricer.jl` | `pricer` for payoffs, arrays, positions, smiles, vectorial processes, market data sets |
| `metrics/pricer_cv.jl` | control variates: `pricer_control_variates`, `control_price` |
| `metrics/distribution.jl` | `distribution` |
| `metrics/delta.jl`, `metrics/rho.jl` | bump-and-reprice with `@set` |
| `metrics/confinterval.jl` | `confinter`, Student intervals |
| `metrics/variance.jl` | `variance`, `independent_samples`, `sample_count`, `per_path_variance`, `payoff_moments` |
| `output_type.jl` | `price_type`, `number_type` |
| `references/metrics.md` | the pipeline, control variates, variance and intervals, output types, known biases |
| `references/configuration.md` | modes, `to_mode`, seeding |

## Process

- Every metric starts with `set_seed!` and moves its inputs once with `to_mode`.
- Results inferred, also for arrays of payoffs of different types (`price_type`, annotated comprehensions).
- Antithetic pairs are one sample in `variance`/`confinter` (#76).

## Outputs

Prices and statistics (one per strike for smiles; `CuArray` in `CudaMode`). Checked by `test/metrics/*.jl`,
`test/methods/control_variates.jl`, `test/models/*.jl`.

## Verify (Kaimon)

Setup: `using FinancialMonteCarlo, Test`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.0))`,
`mc = MonteCarloConfiguration(10_000, 10)`, `opt = EuropeanOption(1.0, 100.0)`.

| Check | Expect |
|:--|:--|
| `@inferred pricer(bs, ZeroRate(0.01), mc, opt)` | `Float64` |
| `pricer(bs, ZeroRate(0.01), mc, opt) == pricer(bs, ZeroRate(0.01), mc, opt)` | `true` (seeded) |
| `typeof(pricer(BlackScholesProcess(0.2f0, Underlying(100f0, 0f0)), ZeroRate(0.01f0), mc, EuropeanOption(1f0, 100f0)))` | `Float32` |
