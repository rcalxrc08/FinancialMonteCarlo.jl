# Monte Carlo configuration

A [`MonteCarloConfiguration`](@ref) holds everything about *how* to simulate: the number of paths
`Nsim`, the number of time steps `Nstep`, the Monte Carlo method and the parallel mode (which also
owns the random number generator and its seed).

```julia
MonteCarloConfiguration(Nsim, Nstep)                                      # standard MC, serial, seed 0
MonteCarloConfiguration(Nsim, Nstep, 42)                                  # seed 42
MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC())  # antithetic variates
MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.MultiThreading())  # parallel mode
MonteCarloConfiguration(Nsim, Nstep, method, parallelMode)                # both
```

## Methods

| Method | Description |
|:--|:--|
| `StandardMC()` | plain Monte Carlo (default) |
| `AntitheticMC()` | antithetic variates: half of the paths use the negated normals of the other half. `Nsim` must be even |
| `SobolMode()` | quasi Monte Carlo with Sobol sequences |
| `ControlVariates(variate, conf, method = StandardMC())` | control variates: a pilot run with configuration `conf` estimates the optimal coefficient and the price of `variate` (exact for a [`Forward`](@ref)), then the main run uses `method`. `variate` must have the same maturity as the priced payoff |
| `PrescribedMC()` | reserved for prescribed random numbers; not implemented by the built-in models |

## Parallel modes and random number generators

| Mode | Description |
|:--|:--|
| `SerialMode(seed = 0, rng = MersenneTwister())` | single thread; any `Random.AbstractRNG` can be used |
| `MultiThreading(sub_mode = SerialMode(), nbatches = Threads.nthreads(), seeds = sub_mode.seed .+ (1:nbatches))` | the paths are split into `nbatches` batches priced on different threads, each with its own seed and its own copy of the generator |
| `MultiProcess(sub_mode = SerialMode(), nbatches = nworkers(), seeds = sub_mode.seed .+ (1:nbatches))` | the same batches, priced on the workers of `Distributed` (needs `using Distributed`, and `@everywhere using FinancialMonteCarlo`) |
| `CudaMode(seed_cpu = 0, seed_gpu = 0)` | simulation on the GPU (needs `using CUDA`) |

Both use exactly `Nsim` paths, and their results depend on the seeds and the number of batches but
not on the number of threads or workers. Each batch is a serial run: prices are the averages of the
batches weighted by their paths, and `variance` and `confinter` are pooled exactly from the
statistics of the batches, so `confinter` is centred on `pricer`. They need a generator that can be
copied, which excludes VectorizedRNG, do not support `SobolMode`, and among the automatic
differentiation backends support only the forward ones.

`MultiProcess` also needs the default `init_rng = true`, because the workers receive copies of the
generators. Functions to differentiate must exist on the workers (`@everywhere f(x) = ...`).
`simulate` and `distribution` send all the paths back to the main process, which can take longer
than simulating them; the other metrics send back only their results.

## Seeds and reproducibility

Every metric resets the random number generator to the seed stored in the configuration before
simulating, so calling a metric twice gives the same result, and bump-and-reprice greeks use
common random numbers. To control the stream yourself, build the generator outside, pass it to
`SerialMode` and either call `simulate` directly or disable the reset with the last argument of
`MonteCarloConfiguration(Nsim, Nstep, method, parallelMode, false)`.

See [Monte Carlo configuration](@ref api-configuration) in the API reference.
