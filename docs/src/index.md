# FinancialMonteCarlo.jl

It currently contains the following capabilities:

- Support for the following Models:
    - Black Scholes
    - Kou
    - Merton
    - Normal Inverse Gaussian
    - Variance Gamma
    - Heston
    - LogNormal Mixture
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
Pkg.clone("https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl.git")
```
## How to Test
After the installation, to test the package type on the Julia REPL the following:
```julia
Pkg.test("FinancialMonteCarlo")
```
## Example of Usage: Pricing Options in Black Scholes Model
The following example shows how to price different kind of options with underlying varying according to the Black Scholes Model, given the implied volatility.
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
spotData1=equitySpotData(S0,r,d);

#Define The Options
Fwd_payoff=Forward(T)
EU_payoff=EuropeanOption(T,K)
AM_payoff=AmericanOption(T,K)
Barrier_payoff=BarrierOptionDownOut(T,K,D)
AsianFloatingStrike_payoff=AsianFloatingStrikeOption(T)
AsianFixedStrike_payoff=AsianFixedStrikeOption(T,K)

#Define the Model
Model=BlackScholesProcess(σ);

#Price
@show FwdPrice=pricer(Model,spotData1,mcConfig,Fwd_payoff);						
@show EuPrice=pricer(Model,spotData1,mcConfig,EU_payoff);
@show AmPrice=pricer(Model,spotData1,mcConfig,AM_payoff);
@show BarrierPrice=pricer(Model,spotData1,mcConfig,Barrier_payoff);
@show AsianPrice1=pricer(Model,spotData1,mcConfig,AsianFloatingStrike_payoff);
@show AsianPrice2=pricer(Model,spotData1,mcConfig,AsianFixedStrike_payoff);
```


## Example of Interaction with DifferentialEquations.jl: Pricing Options in Black Scholes Model
The following example shows how to price different kind of options with underlying varying according to the Black Scholes Model, simulated using DifferentialEquations.jl.
```julia
#Import the Package
using FinancialMonteCarlo,DifferentialEquations;
#Define Spot Datas
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;
u0=S0;
#Define FinancialMonteCarlo Parameters
Nsim=10000;
Nstep=30;
#Define Model Parameters
σ=0.2;
#Build the Structs
mcConfig=MonteCarloConfiguration(Nsim,Nstep);
spotData1=equitySpotData(S0,r,d);

#Define The Options
Fwd_payoff=Forward(T)
EU_payoff=EuropeanOption(T,K)
AM_payoff=AmericanOption(T,K)
Barrier_payoff=BarrierOptionDownOut(T,K,D)
AsianFloatingStrike_payoff=AsianFloatingStrikeOption(T)
AsianFixedStrike_payoff=AsianFixedStrikeOption(T,K)

##Define the Model
#Drift
f(u,p,t) = (r-d)*u
#Diffusion
g(u,p,t) = σ*u
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f,g,u0,tspan)
Model = MonteCarloProblem(prob)
#Price
@show FwdPrice=pricer(Model,spotData1,mcConfig,Fwd_payoff);						
@show EuPrice=pricer(Model,spotData1,mcConfig,EU_payoff);
@show AmPrice=pricer(Model,spotData1,mcConfig,AM_payoff);
@show BarrierPrice=pricer(Model,spotData1,mcConfig,Barrier_payoff);
@show AsianPrice1=pricer(Model,spotData1,mcConfig,AsianFloatingStrike_payoff);
@show AsianPrice2=pricer(Model,spotData1,mcConfig,AsianFixedStrike_payoff);
```
## Keep in mind
There are few things that you should keep in mind when using this library:
- First Order Automatic Differentiation is enabled for any kind of option, also for such that there is no sense (e.g. Binary, Barriers).
- Second Order Automatic Differentiation is enabled for any kind of option but the results are useless.
- Automatic Differentiation is enabled but does not work for process that rely on the simulation of complicated random numbers distributions. If you try you will get a runtime exception.
- Support of Jump Diffusion Differential Equations from DifferentialEquations.jl is broken since [#298](https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/298)


## Functions

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff,mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

## Index

```@index
```