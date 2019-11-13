using Test, FinancialMonteCarlo;
@show "ShiftedLogNormalMixture Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=[0.2,0.2];
lam=Float64[0.999999999];
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8

spotData1=ZeroRateCurve(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=ShiftedLogNormalMixture(sigma,lam,0.0,Underlying(S0,d));

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


@show "Test LogNormalMixture Parameters"
eta=[0.2,0.2]
lam11=[0.9999]
@test_throws(ErrorException,ShiftedLogNormalMixture(-0.2*ones(3),0.1*ones(2),0.0,Underlying(S0,d)))
@test_throws(ErrorException,ShiftedLogNormalMixture(0.2*ones(3),-0.1*ones(2),0.0,Underlying(S0,d)))
@test_throws(ErrorException,ShiftedLogNormalMixture(0.2*ones(3),ones(2),0.0,Underlying(S0,d)))
@test_throws(ErrorException,ShiftedLogNormalMixture(0.2*ones(3),0.2*ones(3),0.0,Underlying(S0,d)))

@test_throws(ErrorException,simulate(ShiftedLogNormalMixture(eta,lam11,0.0,Underlying(S0,d)),spotData1,mc,-T));
