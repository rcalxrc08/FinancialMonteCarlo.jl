using Test, FinancialMonteCarlo;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

underlying_name="ENJ"
underlying_name_2="ABPL"

underlying_=underlying_name*"_"*underlying_name_2

Nsim=10000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8

rfCurve=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model_enj=BlackScholesProcess(sigma,Underlying(S0,d));
Model_abpl=BlackScholesProcess(sigma,Underlying(S0,d));
Model_n=BlackScholesProcess(sigma,Underlying(S0,d));

Model=GaussianCopulaNVariateLogProcess(Model_enj,Model_abpl,0.0)

display(Model)

mktdataset=underlying_|>Model
mktdataset_2=underlying_name|>Model_enj
mktdataset+="notneeded"|>Model_n

portfolio=underlying_name|>1.0*EUData
portfolio+=underlying_name|>1.0*AMData
portfolio+=underlying_name|>-1.0*(-1.0)*BarrierData

price_mkt=pricer(mktdataset,rfCurve,mc,portfolio)
price_old=pricer(mktdataset_2,rfCurve,mc,portfolio)

@test abs(price_mkt-price_old)<1e-8