using Test, FinancialMonteCarlo,DifferentialEquations;
@show "Differential Equation Junctor"
Nsim=10000;
Nstep=30;
toll=0.8;
mc=MonteCarloConfiguration(Nsim,Nstep);
S0=100.0
K=100.0;
D=90.0;
r=0.02
sigma=0.2
T=1.0;
d=0.01;
spotData1=equitySpotData(S0,r,d);
u0=0.0;
#Drift
f1(u,p,t) = (r-d-sigma*sigma/2.0)
#Diffusion
g1(u,p,t) = sigma
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f1,g1,u0,tspan)
monte_prob = MonteCarloProblem(prob)


FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)


@show FwdPrice=pricer(monte_prob,spotData1,mc,FwdData);						
@show EuPrice=pricer(monte_prob,spotData1,mc,EUData);
@show AmPrice=pricer(monte_prob,spotData1,mc,AMData);
@show BarrierPrice=pricer(monte_prob,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(monte_prob,spotData1,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(monte_prob,spotData1,mc,AsianFixedStrikeData);

@test abs(FwdPrice-98.81221412555766)<toll