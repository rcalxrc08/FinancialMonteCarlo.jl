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
ParamDictLongest=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1,"theta1" => theta1, "k1" => k1,"theta11" => theta1, "k11" => k1,"theta12" => theta1, "k12" => k1, "k123" => k1, "k12k" => k1)
ParamDictNeg=Dict{String,Number}("sigma"=>-sigma)
McNeg=MonteCarloBaseData(ParamDictNeg,Nsim,Nstep);
Mc=MonteCarloBaseData(ParamDictLongest,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,McNeg,T));
@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,Mc,-T));
for type1 in concrete_types(AbstractMonteCarloProcess)
	@test_throws(ErrorException,simulate(type1(),spotData1,Mc,T));
end