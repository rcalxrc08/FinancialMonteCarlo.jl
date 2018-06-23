using Base.Test, MonteCarlo,DualNumbers;
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
rho=0.0;
ParamDict=Dict{String,Number}("sigma"=>sigma,"theta"=>theta, "kappa"=>kappa, "sigma_zero" => sigma_zero, "lambda" => lambda, "rho" => rho)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AmericanStdOption(T,K)
BarrierData=BarrierOptionDownOutData(T,K,D)
AsianData1=AsianFloatingStrikeOptionData(T)
AsianData2=AsianFixedStrikeOptionData(T,K)
Model=HestonProcess();

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



@show FwdPrice=pricer(Model,spotData1,mc,FwdData,true,MonteCarlo.antithetic);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,true,MonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,true,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,true,MonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianData1,true,MonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianData2,true,MonteCarlo.antithetic);
tollanti=0.8
@test abs(FwdPrice-98.72567723404445)<tollanti
@test abs(EuPrice-17.62536090688433)<tollanti
@test abs(BarrierPrice-11.38748933756886)<tollanti
@test abs(AsianPrice1-9.762160560168732)<tollanti