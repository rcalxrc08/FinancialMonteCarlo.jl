# FinancialMonteCarlo.jl <img src="etc/logo.png" width="40">  

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rcalxrc08.gitlab.io/FinancialMonteCarlo.jl/)
[![Build Status](https://travis-ci.com/rcalxrc08/FinancialMonteCarlo.jl.svg?branch=master)](https://travis-ci.com/rcalxrc08/FinancialMonteCarlo.jl)[![pipeline status](https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl/badges/master/pipeline.svg)](https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl/commits/master) [![codecov](https://codecov.io/gl/rcalxrc08/FinancialMonteCarlo.jl/branch/master/graph/badge.svg?token=PypZinZKqB)](https://codecov.io/gl/rcalxrc08/FinancialMonteCarlo.jl)
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
    - Shifted LogNormal Mixture
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
mcConfig=MonteCarloConfiguration(Nsim,Nstep); #Configurator
zeroRate=ZeroRate(r);
underlying=Underlying(S0,d); #Underlying relative data

#Define The Option
EU_payoff=EuropeanOption(T,K)
#Define the Model
Model=BlackScholesProcess(σ,underlying);

#Price
@show EuPrice=pricer(Model,zeroRate,mcConfig,EU_payoff);
```

## Curve Support
Non constant zero rates and dividend are managed.
An implied curve is built at time zero, such implied curve is able to return the right implied zero/dividend at a given time,
Without the need to carry the integral structure of the curve.
No support for multicurrency.

## Contracts Algebra
Contracts that refer to the same underlying can be sum togheter in order to build a "new instrument".
In this sense assuming the same underlying, the set of contracts form a vectorial space over "Real" Numbers.

## Market Data and Multivariate Support
A market data set is a dictionary containing all of the process for which we have a view (or a model). ( "KEY" => MODEL)
The portofolio is a dictionary as well but it carries the structure of the portfolio. ( "KEY" => CONTRACTS_ON_MODEL)

## Keep in mind
There are few things that you should keep in mind when using this library:
- First Order Automatic Differentiation is enabled for any kind of option, also for such that there is no sense (e.g. Binary, Barriers).
- Second Order Automatic Differentiation is enabled for any kind of option but the results are useless most of the time.
