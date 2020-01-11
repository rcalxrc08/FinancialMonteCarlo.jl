using Test, FinancialMonteCarlo,DualNumbers;
@show "HestonModel"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
sigma_zero=dual(0.2,1.0);
kappa=0.01;
theta=0.03;
lambda=0.01;
rho_=0.0;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8;

rfCurve=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData1=AsianFloatingStrikeOption(T)
AsianData2=AsianFixedStrikeOption(T,K)
Model=HestonProcess(sigma,sigma_zero,lambda,kappa,rho_,theta,Underlying(S0,d));

@show FwdPrice=pricer(Model,rfCurve,mc,FwdData);
@show EuPrice=pricer(Model,rfCurve,mc,EUData);
@show AmPrice=pricer(Model,rfCurve,mc,AMData);
@show BarrierPrice=pricer(Model,rfCurve,mc,BarrierData);
@show AsianPrice1=pricer(Model,rfCurve,mc,AsianData1);
@show AsianPrice2=pricer(Model,rfCurve,mc,AsianData2);

@test abs(FwdPrice-98.72567723404445)<toll
@test abs(EuPrice-17.62536090688433)<toll
@test abs(BarrierPrice-11.38748933756886)<toll
@test abs(AsianPrice1-9.762160560168732)<toll



@show FwdPrice=pricer(Model,rfCurve,mc1,FwdData);
@show EuPrice=pricer(Model,rfCurve,mc1,EUData);
@show AmPrice=pricer(Model,rfCurve,mc1,AMData);
@show BarrierPrice=pricer(Model,rfCurve,mc1,BarrierData);
@show AsianPrice1=pricer(Model,rfCurve,mc1,AsianData1);
@show AsianPrice2=pricer(Model,rfCurve,mc1,AsianData2);
tollanti=0.8
@test abs(FwdPrice-98.72567723404445)<tollanti
@test abs(EuPrice-17.62536090688433)<tollanti
@test abs(BarrierPrice-11.38748933756886)<tollanti
@test abs(AsianPrice1-9.762160560168732)<tollanti