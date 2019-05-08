using Test, FinancialMonteCarlo,ArrayFire;
ArrayFire.set_seed(UInt64(1)); @show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8

spotData1=equitySpotData(S0,r,d);
FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma);

display(Model)

ArrayFire.set_seed(UInt64(1)); @show FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());

@test abs(FwdPrice-99.1078451563562)<toll
@test abs(EuPrice-8.43005524824866)<toll
@test abs(AmPrice-8.450489415187354)<toll
@test abs(BarrierPrice-7.5008664470880735)<toll
@test abs(AsianPrice1-4.774451704549382)<toll




ArrayFire.set_seed(UInt64(1)); @show FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.antithetic);
ArrayFire.set_seed(UInt64(1)); @show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,FinancialMonteCarlo.antithetic);
ArrayFire.set_seed(UInt64(1)); @show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,FinancialMonteCarlo.antithetic);
ArrayFire.set_seed(UInt64(1)); @show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData,FinancialMonteCarlo.antithetic);
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

ArrayFire.set_seed(UInt64(1)); @show EuPrice=pricer(Model,spotData1,mc,EUDataPut,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AmPrice=pricer(Model,spotData1,mc,AMDataPut,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show BarrierPrice=pricer(Model,spotData1,mc,BarrierDataPut,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeDataPut,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeDataPut,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
tollPut=0.6;
@test abs(FwdPrice-99.1078451563562)<tollPut
@test abs(EuPrice-7.342077422567968)<tollPut
@test abs(AmPrice-7.467065889603365)<tollPut
@test abs(BarrierPrice-0.29071929104261723)<tollPut
@test abs(AsianPrice1-4.230547012372306)<tollPut
@test abs(AsianPrice2-4.236220218194027)<tollPut

#Check Nsim even for antithetic
mc2=MonteCarloConfiguration(Nsim+1,Nstep);
@test_throws(ErrorException,pricer(Model,spotData1,mc2,AsianFixedStrikeData,FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode()));


ArrayFire.set_seed(UInt64(1)); @show var_std_=variance(Model,spotData1,mc,EUDataPut,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode())
ArrayFire.set_seed(UInt64(1)); @show var_ant_=variance(Model,spotData1,mc,EUDataPut,FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode());
@test var_std_>var_ant_
ArrayFire.set_seed(UInt64(1)); @show var_std_=variance(Model,spotData1,mc,[EUDataPut,AMDataPut],FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode())
ArrayFire.set_seed(UInt64(1)); @show var_ant_=variance(Model,spotData1,mc,[EUDataPut,AMDataPut],FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode())
@test prod(var_std_.>=var_ant_);



##############################################################

ArrayFire.set_seed(UInt64(1)); @show IC1=confinter(Model,spotData1,mc,EUDataPut,0.99,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show IC2=confinter(Model,spotData1,mc,EUDataPut,0.99,FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode());
L1=IC1[2]-IC1[1]
L2=IC2[2]-IC2[1]
@test L1>L2

ArrayFire.set_seed(UInt64(1)); @show IC1=confinter(Model,spotData1,mc,[EUDataPut,AMDataPut],0.99,FinancialMonteCarlo.standard,FinancialMonteCarlo.AFMode());
ArrayFire.set_seed(UInt64(1)); @show IC2=confinter(Model,spotData1,mc,[EUDataPut,AMDataPut],0.99,FinancialMonteCarlo.antithetic,FinancialMonteCarlo.AFMode());
L1=IC1[2][2]-IC1[2][1]
L2=IC2[2][2]-IC2[2][1]
@test L1>L2