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
ParamDictBlack=Dict{String,Number}("sigma"=>sigma)
McNeg=MonteCarloBaseData(ParamDictNeg,Nsim,Nstep);
Mc=MonteCarloBaseData(ParamDictLongest,Nsim,Nstep);
McBlack=MonteCarloBaseData(ParamDictBlack,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);


@show "Test Parameters General"
for type1 in concrete_types(AbstractMonteCarloProcess)
	@test_throws(ErrorException,simulate(type1(),spotData1,Mc,T));
end

@show "Test Black Scholes Parameters"
@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,McNeg,T));
@test_throws(ErrorException,simulate(BlackScholesProcess(),spotData1,McBlack,-T));

@show "Test Kou Parameters"
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
ParamDictKou=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => p, "lambdap" => lamp, "lambdam" => lamm)
ParamDictKou1=Dict{String,Number}("sigma"=>sigma, "lambda" => -lam, "p" => p, "lambdap" => lamp, "lambdam" => lamm)
ParamDictKou2=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => -p, "lambdap" => lamp, "lambdam" => lamm)
ParamDictKou3=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => p, "lambdap" => -lamp, "lambdam" => lamm)
ParamDictKou4=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => p, "lambdap" => lamp, "lambdam" => -lamm)
McKou=MonteCarloBaseData(ParamDictKou,Nsim,Nstep);
McKou1=MonteCarloBaseData(ParamDictKou1,Nsim,Nstep);
McKou2=MonteCarloBaseData(ParamDictKou2,Nsim,Nstep);
McKou3=MonteCarloBaseData(ParamDictKou3,Nsim,Nstep);
McKou4=MonteCarloBaseData(ParamDictKou4,Nsim,Nstep);
@test_throws(ErrorException,simulate(KouProcess(),spotData1,McKou,-T));
@test_throws(ErrorException,simulate(KouProcess(),spotData1,McKou1,T));
@test_throws(ErrorException,simulate(KouProcess(),spotData1,McKou2,T));
@test_throws(ErrorException,simulate(KouProcess(),spotData1,McKou3,T));
@test_throws(ErrorException,simulate(KouProcess(),spotData1,McKou4,T));



##### Test Merton
@show "Test Merton Parameters"
mu1=0.03; 
sigma1=0.02;
ParamDictMerton=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "muJ" => mu1, "sigmaJ" => sigma1)
ParamDictMerton1=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "muJ" => mu1, "sigmaJ" => -sigma1)
ParamDictMerton2=Dict{String,Number}("sigma"=>sigma, "lambda" => -lam, "muJ" => mu1, "sigmaJ" => sigma1)
McMerton=MonteCarloBaseData(ParamDictMerton,Nsim,Nstep);
McMerton1=MonteCarloBaseData(ParamDictMerton1,Nsim,Nstep);
McMerton2=MonteCarloBaseData(ParamDictMerton2,Nsim,Nstep);

@test_throws(ErrorException,simulate(MertonProcess(),spotData1,McMerton,-T));
@test_throws(ErrorException,simulate(MertonProcess(),spotData1,McMerton1,T));
@test_throws(ErrorException,simulate(MertonProcess(),spotData1,McMerton2,T));

##### Test Variance Gamma
@show "Test Variance Gamma Parameters"
ParamDictVG=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
ParamDictVG1=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => -k1)
ParamDictVG2=Dict{String,Number}("sigma"=>sigma, "theta" => 1000.0, "k" => k1)
McVG=MonteCarloBaseData(ParamDictVG,Nsim,Nstep);
McVG1=MonteCarloBaseData(ParamDictVG1,Nsim,Nstep);
McVG2=MonteCarloBaseData(ParamDictVG2,Nsim,Nstep);
@test_throws(ErrorException,simulate(VarianceGammaProcess(),spotData1,McVG,-T));
@test_throws(ErrorException,simulate(VarianceGammaProcess(),spotData1,McVG1,T));
@test_throws(ErrorException,simulate(VarianceGammaProcess(),spotData1,McVG2,T));

##### Test NIG
@show "Test NIG Parameters"
ParamDictNIG=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
ParamDictNIG1=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => -k1)
ParamDictNIG2=Dict{String,Number}("sigma"=>sigma, "theta" => 1000.0, "k" => k1)
McNIG=MonteCarloBaseData(ParamDictNIG,Nsim,Nstep);
McNIG1=MonteCarloBaseData(ParamDictNIG1,Nsim,Nstep);
McNIG2=MonteCarloBaseData(ParamDictNIG2,Nsim,Nstep);
@test_throws(ErrorException,simulate(NormalInverseGaussianProcess(),spotData1,McNIG,-T));
@test_throws(ErrorException,simulate(NormalInverseGaussianProcess(),spotData1,McNIG1,T));
@test_throws(ErrorException,simulate(NormalInverseGaussianProcess(),spotData1,McNIG2,T));


##### Test Heston
@show "Test Heston Parameters"
sigma_zero=0.2;
kappa=0.01;
theta=0.03;
lambda=0.01;
rho=0.0;
ParamDictHeston=Dict{String,Number}("sigma"=>sigma,"theta"=>theta, "kappa"=>kappa, "sigma_zero" => sigma_zero, "lambda" => lambda, "rho" => rho)
ParamDictHeston1=Dict{String,Number}("sigma"=>-sigma,"theta"=>theta, "kappa"=>kappa, "sigma_zero" => sigma_zero, "lambda" => lambda, "rho" => rho)
ParamDictHeston2=Dict{String,Number}("sigma"=>sigma,"theta"=>theta, "kappa"=>kappa, "sigma_zero" => -sigma_zero, "lambda" => lambda, "rho" => rho)
ParamDictHeston3=Dict{String,Number}("sigma"=>sigma,"theta"=>theta, "kappa"=>kappa, "sigma_zero" => sigma_zero, "lambda" => lambda, "rho" => -5.0)
ParamDictHeston4=Dict{String,Number}("sigma"=>sigma,"theta"=>theta, "kappa"=>1e-16, "sigma_zero" => sigma_zero, "lambda" => 1e-16, "rho" => rho)
McHeston=MonteCarloBaseData(ParamDictHeston,Nsim,Nstep);
McHeston1=MonteCarloBaseData(ParamDictHeston1,Nsim,Nstep);
McHeston2=MonteCarloBaseData(ParamDictHeston2,Nsim,Nstep);
McHeston3=MonteCarloBaseData(ParamDictHeston3,Nsim,Nstep);
McHeston4=MonteCarloBaseData(ParamDictHeston4,Nsim,Nstep);
@test_throws(ErrorException,simulate(HestonProcess(),spotData1,McHeston,-T));
@test_throws(ErrorException,simulate(HestonProcess(),spotData1,McHeston1,T));
@test_throws(ErrorException,simulate(HestonProcess(),spotData1,McHeston2,T));
@test_throws(ErrorException,simulate(HestonProcess(),spotData1,McHeston3,T));
@test_throws(ErrorException,simulate(HestonProcess(),spotData1,McHeston4,T));



##### Test DifferentialEquations Junctor for Spot Price
toll=0.8;
dict=Dict{String,Number}()
mcDiff=MonteCarloBaseData(dict,Nsim,Nstep);
u0=S0;
#Drift
f(u,p,t) = (r-d)*u
#Diffusion
g(u,p,t) = sigma*u
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f,g,u0,tspan)
monte_prob = MonteCarloProblem(prob)
spotData1=equitySpotData(S0,r,d);

@test_throws(ErrorException,simulate(monte_prob,spotData1,Mc,T));
Tneg=-T;
tspanNeg = (0.0,Tneg)
prob = SDEProblem(f,g,u0,tspanNeg)
monte_probNeg = MonteCarloProblem(prob)
@test_throws(ErrorException,simulate(monte_probNeg,spotData1,mcDiff,Tneg));

##### Test Subordinated BrownianMotion
dt=T/Nstep;
sub=dt*ones(Nsim,Nstep)
ParamDictSub=Dict{String,Number}("drift"=>0.01,"sigma"=>sigma)
McSub=MonteCarloBaseData(ParamDictSub,Nsim,Nstep);
@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(),spotData1,McSub,-T,sub));
@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(),spotData1,Mc,T,sub));

ParamDictBlackNeg=Dict{String,Number}("drift"=>0.01,"sigma"=>-sigma)
McBlackNeg=MonteCarloBaseData(ParamDictBlackNeg,Nsim,Nstep);
@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(),spotData1,McBlackNeg,T,sub));
subw=dt*ones(1,1);
@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(),spotData1,McSub,T,subw));