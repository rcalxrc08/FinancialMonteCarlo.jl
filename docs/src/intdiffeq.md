# Interaction with [DifferentialEquations.jl](http://en.wikipedia.org/wiki/Mixture_model)

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
rfCurve=ZeroRate(S0,r,d);

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
@show FwdPrice=pricer(Model,rfCurve,mcConfig,Fwd_payoff);						
@show EuPrice=pricer(Model,rfCurve,mcConfig,EU_payoff);
@show AmPrice=pricer(Model,rfCurve,mcConfig,AM_payoff);
@show BarrierPrice=pricer(Model,rfCurve,mcConfig,Barrier_payoff);
@show AsianPrice1=pricer(Model,rfCurve,mcConfig,AsianFloatingStrike_payoff);
@show AsianPrice2=pricer(Model,rfCurve,mcConfig,AsianFixedStrike_payoff);
```

## Remarks
Support for Jump Diffusion Processes from DifferentialEquations.jl is currently broken.