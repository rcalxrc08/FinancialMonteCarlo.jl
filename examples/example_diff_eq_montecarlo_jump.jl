using MonteCarlo,DifferentialEquations

Nsim=10000;
Nstep=30;

const dict=Dict{String,Number}()
mc=MonteCarloBaseData(dict,Nsim,Nstep);
S0=100.0
K=100.0;
D=90.0;
r=0.02
sigma=0.2
T=1.0;
d=0.01;
u0=0.0;
sigma=0.2; 
p1=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0; 

#Drift
f(u,p,t) = (r-d-sigma^2/2.0-lam*(p1/(lamp-1)-(1-p1)/(lamm+1)))
#Diffusion
g(u,p,t) = sigma
#Time step
Dt =T/Nstep;
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f,g,u0,tspan)
rate(u,p,t) = (lam*T)
affect!(integrator) = (integrator.u = integrator.u+((rand()<p1)?MonteCarlo.quantile_exp(lamp,rand()):-MonteCarlo.quantile_exp(lamm,rand())))
jump = ConstantRateJump(rate,affect!)
jump_prob = JumpProblem(prob,Direct(),jump)
monte_prob = MonteCarloProblem(jump_prob)


FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData=AsianFloatingStrikeOptionData(T)
spotData1=equitySpotData(S0,r,d);

FwdPrice=pricer(monte_prob,spotData1,mc,FwdData,Forward());
EuPrice=pricer(monte_prob,spotData1,mc,EUData,EuropeanOption());
AMPrice=pricer(monte_prob,spotData1,mc,AMData,AmericanOption());
BarrierPrice=pricer(monte_prob,spotData1,mc,BarrierData,BarrierOptionDownOut());
AsianPrice=pricer(monte_prob,spotData1,mc,AsianData,AsianFloatingStrikeOption());

@show FwdPrice
@show EuPrice
@show AMPrice
@show BarrierPrice
@show AsianPrice