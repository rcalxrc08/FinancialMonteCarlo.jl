# Benchmarks

## Inputs

| Load | Why |
|:--|:--|
| `benchmarks.jl` | `SUITE`: groups `models`, `methods`, `payoffs`, `markets`, `greeks` (10 000 paths × 50 steps) |
| `run.jl` | tunes and runs the groups given as arguments, prints a markdown table (also the job summary) |
| `references/performance.md` | how to compare, CI runner variance |

## Process

- `julia --project=benchmark -t auto benchmark/run.jl [groups]`, `BENCHMARK_SECONDS` per benchmark.
- A new feature with a cost of its own gets a benchmark; compare commits within one machine (runners vary 10-20%).
- `bench_black_ql.jl` and `wip/` are untracked.

## Outputs

The Benchmarks workflow summary; numbers in commit messages for hot-path changes.
