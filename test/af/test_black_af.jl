using Test, FinancialMonteCarlo,ArrayFire;
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
mc=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AFMode());
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC(),FinancialMonteCarlo.AFMode());
toll=0.8

spotData1=ZeroRateCurve(r);
FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

display(Model)

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




@show FwdPrice=pricer(Model,spotData1,mc1,FwdData);
@show EuPrice=pricer(Model,spotData1,mc1,EUData);
@show AmPrice=pricer(Model,spotData1,mc1,AMData,FinancialMonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc1,BarrierData,FinancialMonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc1,AsianFloatingStrikeData,FinancialMonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc1,AsianFixedStrikeData,FinancialMonteCarlo.antithetic);
tollanti=0.6;
@test abs(FwdPrice-99.1078451563562)<tollanti
@test abs(EuPrice-8.43005524824866)<tollanti
@test abs(AmPrice-8.450489415187354)<tollanti
@test abs(BarrierPrice-7.5008664470880735)<tollanti
@test abs(AsianPrice1-4.774451704549382)<tollanti




EUDataPut=EuropeanOption(T,K,false)
AMDataPut=AmericanOption(T,K,false)
BarrierDataPut=BarrierOptionDownOut(T,K,D,false)
AsianFloatingStrikeDataPut=AsianFloatingStrikeOption(T,false)
AsianFixedStrikeDataPut=AsianFixedStrikeOption(T,K,false)

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

#Check Nsim even for antithetic
mc2=MonteCarloConfiguration(Nsim+1,Nstep);
@test_throws(ErrorException,pricer(Model,spotData1,mc2,AsianFixedStrikeData));


@show var_std_=variance(Model,spotData1,mc,EUDataPut)
@show var_ant_=variance(Model,spotData1,mc1,EUDataPut);
@test var_std_>var_ant_
@show var_std_=variance(Model,spotData1,mc,[EUDataPut,AMDataPut])
@show var_ant_=variance(Model,spotData1,mc1,[EUDataPut,AMDataPut])
#@test prod(var_std_.>=var_ant_);



##############################################################

@show IC1=confinter(Model,spotData1,mc,EUDataPut,0.99);
@show IC2=confinter(Model,spotData1,mc1,EUDataPut,0.99);
L1=IC1[2]-IC1[1]
L2=IC2[2]-IC2[1]
@test L1>L2

@show IC1=confinter(Model,spotData1,mc,[EUDataPut,AMDataPut],0.99);
@show IC2=confinter(Model,spotData1,mc1,[EUDataPut,AMDataPut],0.99);
L1=IC1[2][2]-IC1[2][1]
L2=IC2[2][2]-IC2[2][1]
#@test L1>L2