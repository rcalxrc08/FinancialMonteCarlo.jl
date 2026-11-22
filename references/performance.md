# Performance

_The benchmark suite and how to compare._ (Layer 3 reference.)

- `benchmark/benchmarks.jl` defines `SUITE`, a BenchmarkTools group per area (`models`, `methods`, `payoffs`, `markets`,
  `greeks`; 10 000 paths × 50 steps); `benchmark/run.jl` tunes and runs the groups given as arguments (all by default,
  `BENCHMARK_SECONDS` per benchmark, 2 s by default) and prints a markdown table (median, minimum, memory,
  allocations), also written to the job summary on GitHub.
- `Benchmarks.yml` runs it on every push and pull request to master and by hand (`gh workflow run Benchmarks.yml -f
  groups="payoffs"`); shared runners vary by 10-20%, so compare commits within one machine. The whole suite takes about 4
  minutes there (first run, 5eb3e1c: AMD EPYC 7763, 4 threads, European on Black-Scholes 4.1 ms).
- A new feature with a cost of its own gets a benchmark there.
- Measured facts behind the design: AD kernels with rules (Mooncake reverse on Heston was 3-7x the price without them),
  the Zygote rules of #80 (`integral` 483 → 4 μs, `set_parameters!` 0.56-0.67 ms → 11-19 μs), the jump kernel (a thread
  per element, 4x faster than per path), `solve_3x3` against `lu!` (`autodiff.md`, `models.md`, `payoffs.md`).
- Docs: `docs/src/manual/performance.md` and `autodiff.md` tables say the machine they were taken on.
- Zygote: interpolated non-constant arguments, no `$(args)...` splat (FluxML/Zygote.jl#1752, fixed in Zygote 0.7.19:
  still the fair way to time it on older versions).
