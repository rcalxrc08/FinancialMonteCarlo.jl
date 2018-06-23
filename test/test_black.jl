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
AMData=AmericanStdOption(T,K)
BarrierData=BarrierOptionDownOutData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)
Model=BlackScholesProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);						
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData);

@test abs(FwdPrice-99.1078451563562)<toll
@test abs(EuPrice-8.43005524824866)<toll
@test abs(AmPrice-8.450489415187354)<toll
@test abs(BarrierPrice-7.5008664470880735)<toll
@test abs(AsianPrice1-4.774451704549382)<toll


@show FwdPrice=pricer(Model,spotData1,mc,FwdData,MonteCarlo.antithetic);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,MonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,MonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,MonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,MonteCarlo.antithetic);
tollanti=0.6;
@test abs(FwdPrice-99.1078451563562)<tollanti
@test abs(EuPrice-8.43005524824866)<tollanti
@test abs(AmPrice-8.450489415187354)<tollanti
@test abs(BarrierPrice-7.5008664470880735)<tollanti
@test abs(AsianPrice1-4.774451704549382)<tollanti



EUDataPut=EUOptionData(T,K,false)
AMDataPut=AmericanStdOption(T,K,false)
BarrierDataPut=BarrierOptionDownOutData(T,K,D,false)
AsianFloatingStrikeDataPut=AsianFloatingStrikeOptionData(T,false)
AsianFixedStrikeDataPut=AsianFixedStrikeOptionData(T,K,false)
						
@show EuPrice=pricer(Model,spotData1,mc,EUDataPut);
@show AmPrice=pricer(Model,spotData1,mc,AMDataPut);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierDataPut);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeDataPut);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeDataPut);
tollPut=0.6;
@test abs(FwdPrice-99.1078451563562)<tollPut
@test abs(EuPrice-7.342077422567968)<tollPut
@test abs(AmPrice-7.467065889603365)<tollPut
@test abs(BarrierPrice-0.29071929104261723)<tollPut
@test abs(AsianPrice1-4.230547012372306)<tollPut
@test abs(AsianPrice2-4.236220218194027)<tollPut

mc2=MonteCarloBaseData(ParamDict,Nsim+1,Nstep);
@test_throws(ErrorException,pricer(Model,spotData1,mc2,AsianFixedStrikeData,MonteCarlo.antithetic));


############## Complete Options


EUDataBin=BinaryEuropeanOptionData(T,K)
BinAMData=BinaryAmericanOption(T,K)
BarrierDataDI=BarrierOptionDownInData(T,K,D)
BarrierDataUI=BarrierOptionUpInData(T,K,D)
BarrierDataUO=BarrierOptionUpOutData(T,K,D)
doubleBarrierOptionDownOutData=DoubleBarrierOptionData(T,K,K/10.0,1.2*K)

@show BarrierPrice=pricer(Model,spotData1,mc,BarrierDataDI);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierDataUI);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierDataUO);
@show AmBinPrice=pricer(Model,spotData1,mc,BinAMData);
@show EuBinPrice=pricer(Model,spotData1,mc,EUDataBin);
@show doubleBarrier=pricer(Model,spotData1,mc,doubleBarrierOptionDownOutData);