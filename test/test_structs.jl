using MonteCarlo
@show "VarianceGammaProcess"
S0=100.0;
K=100.0;
Kneg=-100.0;
r=0.02;
T=1.0;
Tneg=-1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
theta1=0.01;
k1=0.03;
sigma1=0.02;
ParamDict=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

@test_throws(ErrorException,equitySpotData(-S0,r,d));
#####################################################################
#Payoff Structs
@test_throws(ErrorException,ForwardData(Tneg));
@test_throws(ErrorException,EUOptionData(Tneg,K));
@test_throws(ErrorException,AMOptionData(Tneg,K));
@test_throws(ErrorException,BarrierOptionData(Tneg,K,D));
@test_throws(ErrorException,AsianFloatingStrikeOptionData(Tneg));

@test_throws(ErrorException,EUOptionData(T,Kneg));
@test_throws(ErrorException,AMOptionData(T,Kneg));
@test_throws(ErrorException,BarrierOptionData(T,Kneg,D));
