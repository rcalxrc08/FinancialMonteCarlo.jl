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
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model_enj=BlackScholesProcess(sigma,Underlying(S0,d));
Model_abpl=BlackScholesProcess(sigma,Underlying(S0,d));
Model_tesl=BlackScholesProcess(sigma,Underlying(S0,d));
rho=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
Model=GaussianCopulaNVariateProcess(rho,Model_enj,Model_abpl,Model_tesl)


###############


p=EuropeanOption(1.0,70.0)+Forward(1.0)
p2=3.0*p;
port_=Array{typeof(p)}(undef,7)
port_[1]=p
port_[2]=p2

pricer(Model,rfCurve,mc,port_)