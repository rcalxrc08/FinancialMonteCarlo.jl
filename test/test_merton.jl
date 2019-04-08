using Test, FinancialMonteCarlo;
@show "MertonProcess"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
lam=5.0; 
mu1=0.03; 
sigma1=0.02;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData=AsianFloatingStrikeOption(T)
Model=MertonProcess(sigma,lam,mu1,sigma1);

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice=pricer(Model,spotData1,mc,AsianData);

@test abs(FwdPrice-99.1188767166039)<toll
@test abs(EuPrice-9.084327245917533)<toll
@test abs(BarrierPrice-7.880881290426765)<toll
@test abs(AsianPrice-5.129020349580892)<toll




@show FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.antithetic);
@show EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,FinancialMonteCarlo.antithetic);
@show AsianPrice=pricer(Model,spotData1,mc,AsianData,FinancialMonteCarlo.antithetic);
tollanti=0.6;
@test abs(FwdPrice-99.1188767166039)<tollanti
@test abs(EuPrice-9.084327245917533)<tollanti
@test abs(BarrierPrice-7.880881290426765)<tollanti
@test abs(AsianPrice-5.129020349580892)<tollanti



@show "Test Merton Parameters"
@test_throws(ErrorException,simulate(MertonProcess(sigma,lam,mu1,sigma1),spotData1,mc,-T));
@test_throws(ErrorException,MertonProcess(-sigma,lam,mu1,sigma1))
@test_throws(ErrorException,MertonProcess(sigma,lam,mu1,-sigma1))
@test_throws(ErrorException,MertonProcess(sigma,-lam,mu1,sigma1))
