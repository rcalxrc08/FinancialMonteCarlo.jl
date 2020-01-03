using Test, FinancialMonteCarlo,Statistics;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=[0.00,0.02];
T=1.0;
d=0.01;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8

rfCurve=ZeroRate(r,T);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
EUBin=BinaryEuropeanOption(T,K)
AMData=AmericanOption(T,K)
AmBin=BinaryEuropeanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

display(Model)

@show FwdPrice=pricer(Model,rfCurve,mc,FwdData);
@show EuPrice=pricer(Model,rfCurve,mc,EUData);
@show EuBinPrice=pricer(Model,rfCurve,mc,EUBin);
@show AmPrice=pricer(Model,rfCurve,mc,AMData);
@show AmBinPrice=pricer(Model,rfCurve,mc,AmBin);
@show BarrierPrice=pricer(Model,rfCurve,mc,BarrierData);
@show AsianPrice1=pricer(Model,rfCurve,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(Model,rfCurve,mc,AsianFixedStrikeData);