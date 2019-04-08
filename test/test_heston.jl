using Test, FinancialMonteCarlo;
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
sigma_zero=0.2;
kappa=0.01;
theta=0.03;
lambda=0.01;
rho=0.0;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData1=AsianFloatingStrikeOption(T)
AsianData2=AsianFixedStrikeOption(T,K)
Model=HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta);

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianData1);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianData2);

@test abs(FwdPrice-98.72567723404445)<toll
@test abs(EuPrice-17.62536090688433)<toll
@test abs(BarrierPrice-11.38748933756886)<toll
@test abs(AsianPrice1-9.762160560168732)<toll

@show FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.antithetic);
@show EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,FinancialMonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianData1,FinancialMonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianData2,FinancialMonteCarlo.antithetic);
tollanti=0.8
@test abs(FwdPrice-98.72567723404445)<tollanti
@test abs(EuPrice-17.62536090688433)<tollanti
@test abs(BarrierPrice-11.38748933756886)<tollanti
@test abs(AsianPrice1-9.762160560168732)<tollanti


@test_throws(ErrorException,simulate(HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta),spotData1,McConfig,-T));
@test_throws(ErrorException,HestonProcess(-sigma,sigma_zero,lambda,kappa,rho,theta))
@test_throws(ErrorException,HestonProcess(sigma,-sigma_zero,lambda,kappa,rho,theta))
@test_throws(ErrorException,HestonProcess(sigma,sigma_zero,lambda,kappa,-5.0,theta))
@test_throws(ErrorException,HestonProcess(sigma,sigma_zero,1e-16,1e-16,rho,theta))