# Stage 06: batched modes

Splitting a batched method into serial batches and joining their results; multi-threading.

## Inputs

| Load | Why |
|:--|:--|
| `batches.jl` | `batch_configs`, `batch_sizes`, `batch_results`, `average_on_batches`, `joined_on_batches`, `join_paths`, the batched `pricer`/`simulate` |
| `multi_threading.jl` | `MultiThreading`, its `run_on_batches` |
| `references/batches.md` | how batches are seeded, pooled and typed |
| `references/extensions.md` | `MultiProcess` lives in the Distributed extension |

## Process

- A new batched mode only adds a `run_on_batches` method; its results go through `batch_results` (type stability).
- Results must not depend on the number of threads or workers.

## Outputs

Batched metrics equal in distribution to the serial ones. Checked by `test/methods/parallel.jl`,
`test/methods/multiprocess.jl`.

## Verify (Kaimon)

Setup: `using FinancialMonteCarlo`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.0))`,
`mc = MonteCarloConfiguration(10_000, 10; mode = FinancialMonteCarlo.MultiThreading(4))`.

| Check | Expect |
|:--|:--|
| `pricer(bs, ZeroRate(0.01), mc, EuropeanOption(1.0, 100.0))` | a `Float64`, identical across runs (8.51162738119407 on 2026-11-22) and independent of the thread count |
