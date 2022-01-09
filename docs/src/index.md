# FinancialMonteCarlo.jl

FinancialMonteCarlo.jl contains the following capabilities:

- Support for the following Single Name Models:
    - Black Scholes
    - Kou
    - Merton
    - Normal Inverse Gaussian
    - Variance Gamma
    - Heston
    - LogNormal Mixture
- Support for Sobol and Vectorized RNG.
- Support for Multivariate processes through Gaussian Copula
- Support for non costant zero rates and dividends
- Support for the following Option families:
    - European Options 
    - Barrier Options
    - Asian Options
    - Bermudan Options (Using Longstaff-Schwartz)
    - American Options (Using Longstaff-Schwartz)
- Partial Support for the following Parallelization:
    - CUDA using  [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)
    - Thread based (Native julia)
    - Process based (Native julia)


It also supports the pricing directly from the definition of the Stochastic Differential Equation, using the package [DifferentiatialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

Currently supports [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl), [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) and [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl)
for Automatic Differentiation (where it makes sense).

## How to Install
To install the package simply type on the Julia REPL the following:
```julia
] add FinancialMonteCarlo
```
## How to Test
After the installation, to test the package type on the Julia REPL the following:
```julia
] test FinancialMonteCarlo
```
## Keep in mind
There are few things that you should keep in mind when using this library:
- First Order Automatic Differentiation is enabled for any kind of option, also for such that there is no sense (e.g. Binary, Barriers).
- Second Order Automatic Differentiation is enabled for any kind of option but the results are useless most of the time.

## Index

```@index
```