# FinancialMonteCarlo.jl <img src="etc/logo.png" width="40">  

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rcalxrc08.gitlab.io/FinancialMonteCarlo.jl/)
[![pipeline status](https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl/badges/master/pipeline.svg)](https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl/commits/master) [![coverage report](https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl/badges/master/coverage.svg)](https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl/commits/master)
##### This is a Julia package containing some useful Financial function for Pricing and Risk Management for Equity products.

It currently contains the following capabilities:

- Support for the following Single Name Models:
    - Black Scholes
    - Kou
    - Merton
    - Normal Inverse Gaussian
    - Variance Gamma
    - Heston
    - LogNormal Mixture
- Support for Multivariate process through Gaussian Copula
- Support for the following Option families:
    - European Options 
    - American Options
    - Barrier Options
    - Asian Options
- Partial Support for the following Parallelization:
    - CUDA using  [CuArrays.jl](https://github.com/JuliaGPU/CuArrays.jl)
    - ArrayFire using  [ArrayFire.jl](https://github.com/JuliaComputing/ArrayFire.jl)

It also supports the pricing directly from the definition of the Stochastic Differential Equation, using the package [DifferentiatialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

Currently supports [Dual Numbers](https://github.com/JuliaDiff/DualNumbers.jl), [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [ReverseDiff](https://github.com/JuliaDiff/ReverseDiff.jl)
for Automatic Differentiation (where it makes sense).

Additionally, there has been added support for GPU Parallel execution using [CuArrays](https://github.com/JuliaGPU/CuArrays.jl), you can find examples of this in the benchmark folder.

## How to Install
To install the package simply type on the Julia REPL the following:
```julia
] add https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl
```
## How to Test
After the installation, to test the package type on the Julia REPL the following:
```julia
] test FinancialMonteCarlo
```
## Hello World: Pricing European Call Option in Black Scholes Model
The following example shows how to price a european call option with underlying varying according to the Black Scholes Model, given the volatility:
```julia
#Import the Package
using FinancialMonteCarlo;

#Define Spot Datas
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;
#Define FinancialMonteCarlo Parameters
Nsim=10000;
Nstep=30;
#Define Model Parameters
σ=0.2;
#Build the Structs
mcConfig=MonteCarloConfiguration(Nsim,Nstep);
spotData=ZeroRateCurve(r);
underlying=Underlying(S0,d);

#Define The Options
EU_payoff=EuropeanOption(T,K)
#Define the Model
Model=BlackScholesProcess(σ,underlying);

#Price
@show EuPrice=pricer(Model,spotData,mcConfig,EU_payoff);
```

## Contracts Algebra
TBC

## Market Data and Multivariate Support
TBC

## Keep in mind
There are few things that you should keep in mind when using this library:
- First Order Automatic Differentiation is enabled for any kind of option, also for such that there is no sense (e.g. Binary, Barriers).
- Second Order Automatic Differentiation is enabled for any kind of option but the results are useless most of the time.
- Automatic Differentiation is enabled but does not work for process that rely on the simulation of complicated random numbers distributions. If you try you will get a runtime exception.
