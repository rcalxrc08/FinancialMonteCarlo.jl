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
| `MultiProcess(sub_mode = SerialMode(), nbatches = nworkers(), seeds = 1:nbatches)` | the paths are split into `nbatches` batches priced on the workers of `Distributed` (needs `using Distributed`), each with its own seed |
| `CudaMode(seed_cpu = 0, seed_gpu = 0)` | simulation on the GPU (needs `using CUDA`) |

`MultiProcess` averages the prices of the batches, each with `Nsim ÷ nbatches` paths.

`MultiThreading` uses exactly `Nsim` paths, and its results depend on the seeds and the number of
batches but not on the number of threads. The price is the average of the batch prices weighted by
their paths; the other metrics, such as `variance` and `simulate`, work on the paths of all the
batches, so for American and Bermudan payoffs, whose exercise strategy is estimated from the paths,
`confinter` is not centred exactly on `pricer`. `MultiThreading` needs a generator that can be
copied, which excludes VectorizedRNG, does not support `SobolMode`, and among the automatic
differentiation backends supports only the forward ones.

## Seeds and reproducibility

Every metric resets the random number generator to the seed stored in the configuration before
simulating, so calling a metric twice gives the same result, and bump-and-reprice greeks use
common random numbers. To control the stream yourself, build the generator outside, pass it to
`SerialMode` and either call `simulate` directly or disable the reset with the last argument of
`MonteCarloConfiguration(Nsim, Nstep, method, parallelMode, false)`.

See [Monte Carlo configuration](@ref api-configuration) in the API reference.
