# MonteCarlo
[![pipeline status](https://gitlab.com/rcalxrc08/MonteCarlo.jl/badges/master/pipeline.svg)](https://gitlab.com/rcalxrc08/MonteCarlo.jl/commits/master)
[![coverage report](https://gitlab.com/rcalxrc08/MonteCarlo.jl/badges/master/coverage.svg)](https://gitlab.com/rcalxrc08/MonteCarlo.jl/commits/master)
##### This is a Julia package containing some useful Financial function for Pricing and Risk Management for Equity products.

It currently contains the following capabilities:

- blsprice : Black & Scholes Price for European Options.
- blkprice : Black Price for European Options.
- blsdelta : Black & Scholes Delta sensitivity for European Options.
- blsgamma : Black & Scholes Gamma sensitivity for European Options.
- blstheta : Black & Scholes Theta sensitivity for European Options.
- blsvega  : Black & Scholes Vega sensitivity for European Options.
- blsrho   : Black & Scholes Rho sensitivity for European Options.
- blslambda: Black & Scholes Lambda sensitivity for European Options.
- blspsi   : Black & Scholes Psi sensitivity for European Options.
- blsvanna : Black & Scholes Vanna sensitivity for European Options.
- blsimpv  : Black & Scholes Implied Volatility for European Options (using [Brent Method](http://blog.mmast.net/brent-julia)).
- blkimpv  : Black Implied Volatility for European Options (using [Brent Method](http://blog.mmast.net/brent-julia)).

Currently supports [Dual Numbers](https://github.com/JuliaDiff/DualNumbers.jl) and [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl)
for Automatic Differentiation.

The module depends on [DifferentiatialEquations.jl](https://github.com/JuliaDiff/DualNumbers.jl).

## How to Install
To install the package simply type on the Julia REPL the following:
```Julia
Pkg.clone("https://gitlab.com/rcalxrc08/MonteCarlo.jl.git")
```
## How to Test
After the installation, to test the package type on the Julia REPL the following:
```Julia
Pkg.test("MonteCarlo")
```
## Example of Usage
The following example is the pricing of a Different Kind of Options with underlying varying
according to the Black Scholes Model, given the implied volatility.
```Julia
#Import the Package
using MonteCarlo;
#Define Spot Datas
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;
#Define MonteCarlo Parameters
Nsim=10000;
Nstep=30;
#Define Model Parameters
sigma=0.2;
ParamDict=Dict{String,Number}("sigma"=>sigma)
#Build the Structs
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
spotData1=equitySpotData(S0,r,d);

#Define The Options
FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)

#Define the Model
Model=BlackScholesProcess();

#Price
@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@show AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption());
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownOut());
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption());
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,AsianFixedStrikeOption());
```