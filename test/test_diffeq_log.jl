using Base.Test, MonteCarlo,DifferentialEquations;
@show "Differential Equation Junctor"
Nsim=10000;
Nstep=30;
toll=0.8;
dict=Dict{String,Number}()
mc=MonteCarloBaseData(dict,Nsim,Nstep);
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


FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)


@show FwdPrice=pricer(monte_prob,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(monte_prob,spotData1,mc,EUData,EuropeanOption());
@show AmPrice=pricer(monte_prob,spotData1,mc,AMData,AmericanOption());
@show BarrierPrice=pricer(monte_prob,spotData1,mc,BarrierData,BarrierOptionDownOut());
@show AsianPrice1=pricer(monte_prob,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption());
@show AsianPrice2=pricer(monte_prob,spotData1,mc,AsianFixedStrikeData,AsianFixedStrikeOption());

@test abs(FwdPrice-98.81221412555766)<toll