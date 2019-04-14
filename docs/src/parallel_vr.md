# Variance Reduction and Parallel Mode

The package provides functionalities regarding variance reduction and parallel mode.

## Parallel Mode
The package provides the BaseMode type to represent the type of parallel processe that is going to simulate:

* `SerialMode`: classical serial mode with no parallelization
* `CudaMode`: support for CUDA using CuArrays.jl
* `AFMode`: support for ArrayFire using ArrayFire.jl


## Variance Reduction
The variance reduction is achieved using antithetic variables approach.