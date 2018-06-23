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
spotData1=equitySpotData(S0,r,d);


FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AmericanStdOption(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)


@show FwdPrice=pricer(monte_prob,spotData1,mc,FwdData);						
@show EuPrice=pricer(monte_prob,spotData1,mc,EUData);
@show AmPrice=pricer(monte_prob,spotData1,mc,AMData);
@show BarrierPrice=pricer(monte_prob,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(monte_prob,spotData1,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(monte_prob,spotData1,mc,AsianFixedStrikeData);

@test abs(FwdPrice-98.8436678850961)<toll
@test abs(EuPrice-8.17917706833563)<toll
@test abs(BarrierPrice-7.08419994877601)<toll
@test abs(AsianPrice1-4.70419181812687)<toll