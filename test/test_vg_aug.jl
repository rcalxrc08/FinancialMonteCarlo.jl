using Test, FinancialMonteCarlo,Statistics;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=[0.00,0.02];
T=1.0;
d=FinancialMonteCarlo.Curve([0.00,0.02],T);
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
theta1=0.01; 
k1=0.03;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8;

rfCurve=FinancialMonteCarlo.ZeroRateCurve2(r,T);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
EUBin=BinaryEuropeanOption(T,K)
AMData=AmericanOption(T,K)
AmBin=BinaryEuropeanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=VarianceGammaProcess(sigma,theta1,k1,Underlying(S0,d));

display(Model)
@show FwdPrice=pricer(Model,ZeroRateCurve(r[end]),mc,FwdData);
@show FwdPrice=pricer(Model,rfCurve,mc,FwdData);
@show FwdPrice=pricer(Model,rfCurve,mc1,FwdData);
@show EuPrice=pricer(Model,rfCurve,mc,EUData);
@show EuBinPrice=pricer(Model,rfCurve,mc,EUBin);
@show AmPrice=pricer(Model,rfCurve,mc,AMData);
@show AmBinPrice=pricer(Model,rfCurve,mc,AmBin);
@show BarrierPrice=pricer(Model,rfCurve,mc,BarrierData);
@show AsianPrice1=pricer(Model,rfCurve,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(Model,rfCurve,mc,AsianFixedStrikeData);