using MonteCarlo,DifferentialEquations

Nsim=10000;
Nstep=30;

mc=MonteCarloBaseData(Nsim,Nstep);
S0=100.0
K=100.0;
D=90.0;
r=0.02
sigma=0.2
T=1.0;
d=0.01;
u0=0.0;
#Drift
f(u,p,t) = (r-d-sigma*sigma/2.0)
#Diffusion
g(u,p,t) = sigma
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f,g,u0,tspan)
monte_prob = MonteCarloProblem(prob)

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData=AsianFloatingStrikeOption(T)
spotData1=equitySpotData(S0,r,d);

FwdPrice=pricer(monte_prob,spotData1,mc,FwdData);
EuPrice=pricer(monte_prob,spotData1,mc,EUData);
AMPrice=pricer(monte_prob,spotData1,mc,AMData);
BarrierPrice=pricer(monte_prob,spotData1,mc,BarrierData);
AsianPrice=pricer(monte_prob,spotData1,mc,AsianData);

@show FwdPrice
@show EuPrice
@show AMPrice
@show BarrierPrice
@show AsianPrice