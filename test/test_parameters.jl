using MonteCarlo
@show "Test Parameters"
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
ParamDict=Dict{String,Number}("sigma"=>sigma)
ParamDictLonger=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
ParamDictNeg=Dict{String,Number}("sigma"=>-sigma)
McNeg=MonteCarloBaseData(ParamDictNeg,Nsim,Nstep);
McLonger=MonteCarloBaseData(ParamDictLonger,Nsim,Nstep);
Mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,McLonger,T));
@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,McNeg,T));
@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,Mc,-T));
for type1 in concrete_types(AbstractMonteCarloProcess)
	if type1()!=BlackScholesProcess()&&type1()!=SubordinatedBrownianMotion()
		@test_throws(ErrorException,simulate(type1(),spotData1,Mc,T));
	end
end