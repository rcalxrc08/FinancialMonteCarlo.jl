using Base.Test, MonteCarlo,DifferentialEquations;

Nsim=10000;
Nstep=30;

dict=Dict{String,Number}()
mc=MonteCarloBaseData(dict,Nsim,Nstep);
S0=100.0
K=100.0;
D=90.0;
r=0.02
sigma=0.2
T=1.0;
d=0.01;
u0=S0;
#Drift
f(u,p,t) = (r-d)*u
#Diffusion
g(u,p,t) = sigma*u
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f,g,u0,tspan)
monte_prob = MonteCarloProblem(prob)

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AmericanStdOption(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData=AsianFloatingStrikeOptionData(T)
spotData1=equitySpotData(S0,r,d);

optionDatas=[FwdData,EUData,AMData,BarrierData,AsianData]
options=[Forward()]

optPrices=pricer(monte_prob,spotData1,mc,optionDatas,options);

@test 0==0