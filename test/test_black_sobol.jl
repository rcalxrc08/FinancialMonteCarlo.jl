using Test, FinancialMonteCarlo;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
d_=FinancialMonteCarlo.Curve([0.009999,0.01],T);
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.SobolMode());
toll=0.8

spotData1=ZeroRate(r);

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

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show EuBinPrice=pricer(Model,spotData1,mc,EUBin);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show AmBinPrice=pricer(Model,spotData1,mc,AmBin);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData);

@test abs(FwdPrice-99.1078451563562)<toll
@test abs(FwdPrice_-99.1078451563562)<toll
@test abs(EuPrice-8.43005524824866)<toll
@test abs(AmPrice-8.450489415187354)<toll
@test abs(BarrierPrice-7.5008664470880735)<toll
@test abs(AsianPrice1-4.774451704549382)<toll