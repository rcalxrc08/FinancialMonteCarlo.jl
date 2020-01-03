using Test, FinancialMonteCarlo;
@show "HestonModel"
S0=100.0;
K=100.0;
r=[0.02;0.02];
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
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8;

rfCurve=ZeroRate(r,T);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData1=AsianFloatingStrikeOption(T)
AsianData2=AsianFixedStrikeOption(T,K)
Model=HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta,Underlying(S0,d));

@show FwdPrice=pricer(Model,rfCurve,mc,FwdData);

@test abs(FwdPrice-98.72567723404445)<toll

@show FwdPrice=pricer(Model,rfCurve,mc1,FwdData);
tollanti=0.8
@test abs(FwdPrice-98.72567723404445)<tollanti


@test_throws(ErrorException,simulate(HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta,Underlying(S0,d)),rfCurve,mc,-T));
@test_throws(ErrorException,simulate(HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta,Underlying(S0,d)),rfCurve,mc1,-T));
@test_throws(ErrorException,HestonProcess(-sigma,sigma_zero,lambda,kappa,rho,theta,Underlying(S0,d)))
@test_throws(ErrorException,HestonProcess(sigma,-sigma_zero,lambda,kappa,rho,theta,Underlying(S0,d)))
@test_throws(ErrorException,HestonProcess(sigma,sigma_zero,lambda,kappa,-5.0,theta,Underlying(S0,d)))
@test_throws(ErrorException,HestonProcess(sigma,sigma_zero,1e-16,1e-16,rho,theta,Underlying(S0,d)))