using Base.Test, MonteCarlo;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2;
ParamDict=Dict{String,Number}("sigma"=>sigma)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=1e-3;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)
Model=BlackScholesProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@show AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption());
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownOut());
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption());
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,AsianFixedStrikeOption());

@test abs(FwdPrice-99.1078451563562)<toll
@test abs(EuPrice-8.43005524824866)<toll
@test abs(AmPrice-8.450489415187354)<toll
@test abs(BarrierPrice-7.5008664470880735)<toll
@test abs(AsianPrice1-4.774451704549382)<toll


@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward(),true,MonteCarlo.antithetic);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption(),true,MonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption(),true,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownOut(),true,MonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption(),true,MonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,AsianFixedStrikeOption(),true,MonteCarlo.antithetic);
tollanti=0.6;
@test abs(FwdPrice-99.1078451563562)<tollanti
@test abs(EuPrice-8.43005524824866)<tollanti
@test abs(AmPrice-8.450489415187354)<tollanti
@test abs(BarrierPrice-7.5008664470880735)<tollanti
@test abs(AsianPrice1-4.774451704549382)<tollanti



@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward(),false);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption(),false);
@show AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption(),false);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownOut(),false);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption(),false);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,AsianFixedStrikeOption(),false);
tollPut=0.6;
@test abs(FwdPrice-99.1078451563562)<tollPut
@test abs(EuPrice-7.342077422567968)<tollPut
@test abs(AmPrice-7.467065889603365)<tollPut
@test abs(BarrierPrice-0.29071929104261723)<tollPut
@test abs(AsianPrice1-4.230547012372306)<tollPut
@test abs(AsianPrice2-4.236220218194027)<tollPut



@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward(),false,MonteCarlo.antithetic);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption(),false,MonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption(),false,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownOut(),false,MonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption(),false,MonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,AsianFixedStrikeOption(),false,MonteCarlo.antithetic);
tollAntiPut=0.6;
@test abs(FwdPrice-99.1078451563562)<tollAntiPut
@test abs(EuPrice-7.383764847221101)<tollAntiPut
@test abs(AmPrice-7.423390795951484)<tollAntiPut
@test abs(BarrierPrice-0.300813355608894)<tollAntiPut
@test abs(AsianPrice1-4.264627139486553)<tollAntiPut
@test abs(AsianPrice2-4.192045764856742)<tollAntiPut


############## Complete Options


FwdData=ForwardData(T)
EUDataBin=BinaryEuropeanOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)
doubleBarrierOptionData=DoubleBarrierOptionData(T,K,K/10.0,1.2*K)
Model=BlackScholesProcess();

@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownIn(),false,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionUpIn(),false,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionUpOut(),false,MonteCarlo.antithetic);
@show AmBinPrice=pricer(Model,spotData1,mc,AMData,BinaryAmericanOption(),false,MonteCarlo.antithetic);
@show EuBinPrice=pricer(Model,spotData1,mc,EUDataBin,BinaryEuropeanOption(),false,MonteCarlo.antithetic);
@show doubleBarrier=pricer(Model,spotData1,mc,doubleBarrierOptionData,DoubleBarrierOption(),false,MonteCarlo.antithetic);