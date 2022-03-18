# Variance Reduction, Parallel Mode and Random Number Generators

The package provides functionalities regarding variance reduction and parallel mode.

## Parallel Mode
The package provides the BaseMode type to represent the type of parallel processe that is going to simulate:

* `SerialMode`: classical serial mode with no parallelization
* `CudaMode`: support for CUDA using CUDA.jl
* `MultiThreading`: Threads parallelization with julia built in costructs
* `MultiProcess`: Processes parallelization with julia built in costructs

## Random Number Generator and seed control
You can provide your own rng (that must inherit from Base.AbstractRNG), or choose the built in one controlling the seed.
Seed is controlled from inside the package, each call to the constructor of MonteCarloConfiguration or pricer(...) automatically set the seed.
If you want to control the seed from outside, just provide already build rng to constructor of MonteCarloConfiguration and avoid calls to pricer(...).

## Variance Reduction
The variance reduction is achieved using antithetic variables approach and control variates.

## How to

In order to impose a particular parallel mode, rng or variance reduction tecnique you should provide additional parameters to MonteCarloConfiguration constructor.