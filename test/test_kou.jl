using Base.Test, MonteCarlo;
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
ParamDict=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => p, "lambdap" => lamp, "lambdam" => lamm)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AmericanStdOption(T,K)
BarrierData=BarrierOptionDownOutData(T,K,D)
AsianData1=AsianFloatingStrikeOptionData(T)
AsianData2=AsianFixedStrikeOptionData(T,K)
Model=KouProcess();

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


@show FwdPrice=pricer(Model,spotData1,mc,FwdData,true,MonteCarlo.antithetic);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,true,MonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,true,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,true,MonteCarlo.antithetic);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianData1,true,MonteCarlo.antithetic);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianData2,true,MonteCarlo.antithetic);
tollanti=0.6
@test abs(FwdPrice-99.41332633109904)<tollanti
@test abs(EuPrice-10.347332240535199)<tollanti
@test abs(BarrierPrice-8.860123655599818)<tollanti
@test abs(AsianPrice1-5.81798437145069)<tollanti