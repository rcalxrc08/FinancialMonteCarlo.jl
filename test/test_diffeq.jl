using Test, FinancialMonteCarlo,DifferentialEquations;
@show "Differential Equation Junctor"
Nsim=10000;
Nstep=30;
toll=0.8;
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
dr_(u,p,t) = (r-d)*u
#Diffusion
g_1(u,p,t) = sigma*u
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(dr_,g_1,u0,tspan)
monte_prob = MonteCarloProblem(prob)
rfCurve=ZeroRateCurve(r);


FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)


@show FwdPrice=pricer(monte_prob,rfCurve,mc,FwdData);
@show EuPrice=pricer(monte_prob,rfCurve,mc,EUData);
@show AmPrice=pricer(monte_prob,rfCurve,mc,AMData);
@show BarrierPrice=pricer(monte_prob,rfCurve,mc,BarrierData);
@show AsianPrice1=pricer(monte_prob,rfCurve,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(monte_prob,rfCurve,mc,AsianFixedStrikeData);

@test abs(FwdPrice-98.8436678850961)<toll
@test abs(EuPrice-8.17917706833563)<toll
@test abs(BarrierPrice-7.08419994877601)<toll
@test abs(AsianPrice1-4.70419181812687)<toll



Tneg=-T;
tspanNeg = (0.0,Tneg)
probNeg = SDEProblem(f,g,u0,tspanNeg)
monte_probNeg = MonteCarloProblem(probNeg)
@test_throws(ErrorException,simulate(monte_probNeg,rfCurve,mc,Tneg));
