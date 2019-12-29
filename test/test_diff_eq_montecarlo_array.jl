using Test, FinancialMonteCarlo,DifferentialEquations;

Nsim=10000;
Nstep=30;

mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
S0=100.0
K=100.0;
D=90.0;
r=0.02
sigma=0.2
T=1.0;
d=0.01;
u0=S0;
#Drift
drift_(u,p,t) = (r-d)*u
#Diffusion
diffusion_(u,p,t) = sigma*u
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(drift_,diffusion_,u0,tspan)
monte_prob = MonteCarloProblem(prob)

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData=AsianFloatingStrikeOption(T)
rfCurve=ZeroRateCurve(r);

optionDatas=[FwdData,EUData,AMData,BarrierData,AsianData]

optPrices=pricer(monte_prob,rfCurve,mc,optionDatas);
@show optPrices
@test 0==0