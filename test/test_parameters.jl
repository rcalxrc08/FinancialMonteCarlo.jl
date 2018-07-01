using MonteCarlo, Base.Test,DifferentialEquations
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
McConfig=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

@show "Test Black Scholes Parameters"
@test_throws(ErrorException,simulate(BlackScholesProcess(sigma),spotData1,McConfig,-T));
@test_throws(ErrorException,BlackScholesProcess(-sigma))

@show "Test LogNormalMixture Parameters"
eta=[0.2,0.2]
lam11=Float64[0.9999]
@test_throws(ErrorException,LogNormalMixture(-0.2*ones(3),0.1*ones(2)))
@test_throws(ErrorException,LogNormalMixture(0.2*ones(3),-0.1*ones(2)))
@test_throws(ErrorException,LogNormalMixture(0.2*ones(3),ones(2)))
@test_throws(ErrorException,LogNormalMixture(0.2*ones(3),0.2*ones(3)))

@test_throws(ErrorException,simulate(LogNormalMixture(eta,lam11),spotData1,McConfig,-T));

@show "Test Kou Parameters"
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;

@test_throws(ErrorException,simulate(KouProcess(sigma, lam, p,  lamp,  lamm),spotData1,McConfig,-T));
@test_throws(ErrorException,KouProcess(-sigma,lam, p,  lamp,  lamm))
@test_throws(ErrorException,KouProcess( sigma, -lam, p,  lamp,  lamm))
@test_throws(ErrorException,KouProcess( sigma, lam,-p,  lamp,  lamm))
@test_throws(ErrorException,KouProcess( sigma, lam, p, -lamp,  lamm))
@test_throws(ErrorException,KouProcess( sigma, lam, p,  lamp, -lamm))



##### Test Merton
@show "Test Merton Parameters"
mu1=0.03; 
sigma1=0.02;


@test_throws(ErrorException,simulate(MertonProcess(sigma,lam,mu1,sigma1),spotData1,McConfig,-T));
@test_throws(ErrorException,MertonProcess(-sigma,lam,mu1,sigma1))
@test_throws(ErrorException,MertonProcess(sigma,lam,mu1,-sigma1))
@test_throws(ErrorException,MertonProcess(sigma,-lam,mu1,sigma1))

##### Test Variance Gamma
@show "Test Variance Gamma Parameters"

k=0.01;
@test_throws(ErrorException,simulate(VarianceGammaProcess(sigma,theta1,k),spotData1,McConfig,-T));
@test_throws(ErrorException,VarianceGammaProcess(-sigma,theta1,k))
@test_throws(ErrorException,VarianceGammaProcess(sigma,theta1,-k))
@test_throws(ErrorException,VarianceGammaProcess(sigma,10000.0,k))

##### Test NIG
@show "Test NIG Parameters"

@test_throws(ErrorException,simulate(NormalInverseGaussianProcess(sigma,theta1,k1),spotData1,McConfig,-T));
@test_throws(ErrorException,NormalInverseGaussianProcess(-sigma,theta1,k1))
@test_throws(ErrorException,NormalInverseGaussianProcess(sigma,theta1,-k1))
@test_throws(ErrorException,NormalInverseGaussianProcess(sigma,10000.0,k1))


##### Test Heston
@show "Test Heston Parameters"
sigma_zero=0.2;
kappa=0.01;
theta=0.03;
lambda=0.01;
rho=0.0;


@test_throws(ErrorException,simulate(HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta),spotData1,McConfig,-T));
@test_throws(ErrorException,HestonProcess(-sigma,sigma_zero,lambda,kappa,rho,theta))
@test_throws(ErrorException,HestonProcess(sigma,-sigma_zero,lambda,kappa,rho,theta))
@test_throws(ErrorException,HestonProcess(sigma,sigma_zero,lambda,kappa,-5.0,theta))
@test_throws(ErrorException,HestonProcess(sigma,sigma_zero,1e-16,1e-16,rho,theta))



##### Test DifferentialEquations Junctor for Spot Price
toll=0.8;

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

Tneg=-T;
tspanNeg = (0.0,Tneg)
prob = SDEProblem(f,g,u0,tspanNeg)
monte_probNeg = MonteCarloProblem(prob)
@test_throws(ErrorException,simulate(monte_probNeg,spotData1,McConfig,Tneg));

##### Test Subordinated BrownianMotion
drift=0.0;
dt=T/Nstep;
sub=dt*ones(Nsim,Nstep)

@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(sigma,drift),spotData1,McConfig,-T,sub));
@test_throws(ErrorException,SubordinatedBrownianMotion(-sigma,drift))


@test_throws(ErrorException,SubordinatedBrownianMotion(-sigma,drift))
subw=dt*ones(1,1);
@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(sigma,drift),spotData1,McConfig,T,subw));