using FinancialMonteCarlo

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=300;
sigma=0.2;
variate_=FinancialMonteCarlo.ControlVariates(Forward(T),MonteCarloConfiguration(1000,100))
variate_2=FinancialMonteCarlo.ControlVariates(AsianFloatingStrikeOption(T),MonteCarloConfiguration(1000,100))
mc=MonteCarloConfiguration(Nsim,Nstep,variate_);
mc2=MonteCarloConfiguration(Nsim,Nstep,variate_);
mc1=MonteCarloConfiguration(Nsim,Nstep);
toll=1e-3;

spotData1=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData);
AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData);
AsianPrice2=pricer(Model,spotData1,mc2,AsianFixedStrikeData);