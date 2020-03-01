using Test, FinancialMonteCarlo,DualNumbers,BenchmarkTools;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

underlying_name="ENJ"
underlying_name_2="ABPL"

underlying_=underlying_name*"_"*underlying_name_2*"_TESL"

Nsim=10000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8

rfCurve=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOptionND(T,K)
Model_enj=BlackScholesProcess(sigma,Underlying(S0,d));
Model_enj_2=BlackScholesProcess(dual(sigma,1.0),Underlying(S0,d));
Model_abpl=BlackScholesProcess(sigma,Underlying(S0,d));
Model_tesl=BlackScholesProcess(sigma,Underlying(S0,d));
rho_1=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
Model=GaussianCopulaNVariateProcess(rho_1,Model_enj,Model_abpl,Model_tesl)
ModelDual=GaussianCopulaNVariateProcess(rho_1,Model_enj_2,Model_abpl,Model_tesl)

display(Model)

mktdataset=underlying_|>Model
mktdataset_dual=underlying_|>ModelDual

portfolio_=[EUData];
portfolio=underlying_|>EUData

@btime price_mkt=pricer(mktdataset,rfCurve,mc,portfolio)
@btime price_mkt=pricer(mktdataset_dual,rfCurve,mc,portfolio)
