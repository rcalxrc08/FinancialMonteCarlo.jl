# Variance Reduction and Parallel Mode

The package provides functionalities regarding variance reduction and parallel mode.

## Parallel Mode
The package provides the BaseMode type to represent the type of parallel processe that is going to simulate:

* `SerialMode`: classical serial mode with no parallelization
* `CudaMode`: support for CUDA using CuArrays.jl
* `AFMode`: support for ArrayFire using ArrayFire.jl


## Variance Reduction
The variance reduction is achieved using antithetic variables approach.



## Parallelization vs Number Type
Keep in mind the following table when trying to use parallelization and automatic differentiation / complex numbers

**Element Type** | **Parallelization mode** | **Descriptions**
--- | --- | ---
`Float X` | `Serial` | Works with any type of Float
`Float X` | `GPU` | Works with any type of Float, best with Float32
`num<:Number` | `Serial` | Works almost everywhere with any type of Number, like Duals,Complex,Hypers,etc.
`num<:Number` | `GPU` | Works just with Complex and ForwardDiff.Dual.
