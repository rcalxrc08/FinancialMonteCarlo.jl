using Base.Test, MonteCarlo;
@show "VarianceGammaProcess"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
theta1=0.01; 
k1=0.03;
ParamDict=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AmericanStdOption(T,K)
BarrierData=BarrierOptionDownOutData(T,K,D)
AsianData=AsianFloatingStrikeOptionData(T)
Model=VarianceGammaProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);						
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice=pricer(Model,spotData1,mc,AsianData);

@test abs(FwdPrice-97.94751460264095)<toll
@test abs(EuPrice-7.738298817933206)<toll
@test abs(BarrierPrice-6.886023820038332)<toll
@test abs(AsianPrice-4.414948846776423)<toll



@show FwdPrice=pricer(Model,spotData1,mc,FwdData,MonteCarlo.antithetic);						
@show EuPrice=pricer(Model,spotData1,mc,EUData,MonteCarlo.antithetic);
@show AmPrice=pricer(Model,spotData1,mc,AMData,MonteCarlo.antithetic);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,MonteCarlo.antithetic);
@show AsianPrice=pricer(Model,spotData1,mc,AsianData,MonteCarlo.antithetic);

@test abs(FwdPrice-97.94751460264095)<toll
@test abs(EuPrice-7.738298817933206)<toll
@test abs(BarrierPrice-6.886023820038332)<toll
@test abs(AsianPrice-4.414948846776423)<toll