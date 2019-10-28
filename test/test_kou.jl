using Test, FinancialMonteCarlo;
@show "KouModel"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8;

spotData1=equitySpotData(r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianData1=AsianFloatingStrikeOption(T)
AsianData2=AsianFixedStrikeOption(T,K)
Model=KouProcess(sigma,lam,p,lamp,lamm,S0);

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianData1);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianData2);

@test abs(FwdPrice-99.41332633109904)<toll
@test abs(EuPrice-10.347332240535199)<toll
@test abs(BarrierPrice-8.860123655599818)<toll
@test abs(AsianPrice1-5.81798437145069)<toll


@show FwdPrice=pricer(Model,spotData1,mc1,FwdData);
@show EuPrice=pricer(Model,spotData1,mc1,EUData);
@show AmPrice=pricer(Model,spotData1,mc1,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc1,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc1,AsianData1);
@show AsianPrice2=pricer(Model,spotData1,mc1,AsianData2);
tollanti=0.6
@test abs(FwdPrice-99.41332633109904)<tollanti
@test abs(EuPrice-10.347332240535199)<tollanti
@test abs(BarrierPrice-8.860123655599818)<tollanti
@test abs(AsianPrice1-5.81798437145069)<tollanti




@show "Test Kou Parameters"

@test_throws(ErrorException,simulate(KouProcess(sigma, lam, p,  lamp,  lamm,S0),spotData1,mc,-T));
@test_throws(ErrorException,KouProcess(-sigma,lam, p,  lamp,  lamm,S0))
@test_throws(ErrorException,KouProcess( sigma, -lam, p,  lamp,  lamm,S0))
@test_throws(ErrorException,KouProcess( sigma, lam,-p,  lamp,  lamm,S0))
@test_throws(ErrorException,KouProcess( sigma, lam, p, -lamp,  lamm,S0))
@test_throws(ErrorException,KouProcess( sigma, lam, p,  lamp, -lamm,S0))
