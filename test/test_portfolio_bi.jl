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

spotData1=ZeroRateCurve(r);

FwdData=Forward(T)
EUData=EuropeanOption2D(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model_enj=BlackScholesProcess(sigma,Underlying(S0,d));
Model_abpl=BlackScholesProcess(sigma,Underlying(S0,d));

Model=GaussianCopulaBivariateProcess(Model_enj,Model_enj,0.0,)

display(Model)

mktdataset=underlying_|>Model

portfolio_=[EUData];
portfolio=underlying_|>EUData

price_mkt=pricer(mktdataset,spotData1,mc,portfolio)
price_old= sum(pricer(Model,spotData1,mc,portfolio_))


@test abs(price_mkt-price_old)<1e-8