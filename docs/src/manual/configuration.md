# Monte Carlo configuration

A [`MonteCarloConfiguration`](@ref) is the Monte Carlo pricing method: the method argument of
[`pricer`](@ref) and of the other metrics. It holds the number of paths `Nsim`, the number of time
steps `Nstep`, and, as keywords, the sampling scheme, the execution mode, the seed and the random
number generator.

```julia
MonteCarloConfiguration(Nsim, Nstep)                                                 # standard MC, serial, seed 0
MonteCarloConfiguration(Nsim, Nstep; seed = 42)                                      # seed 42
MonteCarloConfiguration(Nsim, Nstep; sampling = FinancialMonteCarlo.AntitheticMC())  # antithetic variates
MonteCarloConfiguration(Nsim, Nstep; mode = FinancialMonteCarlo.MultiThreading())    # on all the threads
MonteCarloConfiguration(Nsim, Nstep; rng = Xoshiro(), reseed = false)                # your generator, its stream continues
```

Every combination is checked when the configuration is built, so an unsupported one (antithetic
variates with an odd `Nsim`, Sobol sequences in batches or on the GPU, `MultiProcess` with
`reseed = false`, `CudaMode` without CUDA.jl) fails there and not when pricing.

## Sampling

| Sampling | Description |
|:--|:--|
| `StandardMC()` | plain Monte Carlo (default) |
| `AntitheticMC()` | antithetic variates: half of the paths use the negated normals of the other half (and, for Kou and Merton, the uniforms `1 - u` of its jumps). `Nsim` must be even |
| `SobolMC()` | quasi Monte Carlo with Sobol sequences, only in `SerialMode`: each path is one point of dimension `Nstep` (`2 Nstep` for Heston, `Nstep + 1 + 2M` for Kou and Merton), after skipping the first `2^m - 1` points, `2^m ≤ Nstep · Nsim + 1 < 2^(m+1)` |

With `AntitheticMC` and `SobolMC` the jumps of Kou and Merton are drawn from uniforms by inversion:
one uniform gives the number of jumps (Poisson), then each jump takes one for its time and one for
its size. The antithetic path uses `1 - u` for all of them. With `SobolMC` a path has room for `M`
jumps, the Poisson quantile at `1 - sqrt(eps)`: the rare path with more draws them from `rng`.
| `ControlVariates(variate, conf, sampling = StandardMC())` | control variates: a pilot run with configuration `conf` estimates the optimal coefficient and the price of `variate` (exact for a [`Forward`](@ref)), then the main run uses `sampling`. `variate` must have the same maturity as the priced payoff |

## Execution modes

The mode says where and how the method runs. It holds no seed and no generator, so it is the same
object for the methods without randomness of other packages (the Fourier methods of FinancialFFT,
for example): every pricing method carries its mode, read with `FinancialMonteCarlo.execution_mode(method)`
and changed with `FinancialMonteCarlo.with_mode(method, mode)`.

| Mode | Description |
|:--|:--|
| `SerialMode()` | one thread of the CPU (default) |
| `MultiThreading(nbatches = Threads.nthreads())` | the paths are split into `nbatches` batches priced on different threads |
| `MultiProcess(nbatches = nworkers())` | the same batches, priced on the workers of `Distributed` (needs `using Distributed`, and `@everywhere using FinancialMonteCarlo`) |
| `CudaMode()` | simulation on the GPU (needs `using CUDA`) |

The batched modes use exactly `Nsim` paths, and their results depend on the seeds and the number of
batches but not on the number of threads or workers. Batch `i` has its own copy of the generator,
with seed `seed + i`, or `seed[i]` when `seed` is a vector with one seed per batch. Each batch is a
serial run: prices are the averages of the batches weighted by their paths, and `variance` and
`confinter` are pooled exactly from the statistics of the batches, so `confinter` is centred on
`pricer`. They need a generator that can be copied, which excludes VectorizedRNG, do not support
`SobolMC`, and among the automatic differentiation backends support only the forward ones.

`MultiProcess` also needs the default `reseed = true`, because the workers receive copies of the
generators. Functions to differentiate must exist on the workers (`@everywhere f(x) = ...`).
`simulate` and `distribution` send all the paths back to the main process, which can take longer
than simulating them; the other metrics send back only their results.

In `CudaMode` the paths are simulated on the GPU, from the CURAND generator seeded with `seed`
(the jump times and the subordinators are drawn on the host from `rng`). The mode decides where the
computation runs: a smile with its strikes on the host is moved to the GPU once per metric, and one
with its strikes on the GPU is moved to the host in the other modes. The results stay on the device
that computed them.

## Seeds and reproducibility

Every metric resets the random number generators to the seeds stored in the configuration before
simulating, so calling a metric twice gives the same result, and bump-and-reprice greeks use
common random numbers. To control the stream yourself, pass your generator with `rng` and either
call `simulate` directly or disable the reset with `reseed = false`.

See [Monte Carlo configuration](@ref api-configuration) in the API reference.
