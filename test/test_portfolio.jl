using Test, FinancialMonteCarlo;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

underlying_name="ENJ"

Nsim=10000;
Nstep=30;
sigma=0.2;
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
Model=BlackScholesProcess(sigma,Underlying(S0,d));

display(Model)

mktdataset=underlying_name|>Model
portfolio_=[FwdData;EUData;AMData;BarrierData;AsianFloatingStrikeData;AsianFixedStrikeData];
portfolio=underlying_name|>FwdData
portfolio+=underlying_name|>EUData
portfolio+=underlying_name|>AMData
portfolio+=underlying_name|>BarrierData
portfolio+=underlying_name|>AsianFloatingStrikeData
portfolio+=underlying_name|>AsianFixedStrikeData


price_mkt=pricer(mktdataset,spotData1,mc,portfolio)
price_old= sum(pricer(Model,spotData1,mc,portfolio_))


@test abs(price_mkt-price_old)<1e-8