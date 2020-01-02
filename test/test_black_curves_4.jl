using Test, FinancialMonteCarlo,Statistics;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=[0.00,0.02];
T=1.0;
d=FinancialMonteCarlo.Curve([0.019,0.02],T);
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=1.8

rfCurve=FinancialMonteCarlo.ZeroRateCurve2(r,T);

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

@test abs(FwdPrice-99.1078451563562)<toll
@test abs(EuPrice-8.43005524824866)<toll
@test abs(AmPrice-8.450489415187354)<toll
@test abs(BarrierPrice-7.5008664470880735)<toll
@test abs(AsianPrice1-4.774451704549382)<toll