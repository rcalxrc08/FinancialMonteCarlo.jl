using Test, FinancialMonteCarlo,DifferentialEquations;
quantile_exp(lam::Number,rand1::Number)::Number=-log(1-rand1)/lam;
@show "Differential Equation Junctor"
Nsim=10000;
Nstep=30;
sigma=0.2; 
p1=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
toll=0.8
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
S0=100.0
K=100.0;
D=90.0;
r=0.02
sigma=0.2
T=1.0;
d=0.01;
spotData1=equitySpotData(S0,r,d);
u0=0.0;
#Drift
f1(u,p,t) = (r-d-sigma^2/2.0-lam*(p1/(lamp-1)-(1-p1)/(lamm+1)))
#Diffusion
g1_(u,p,t) = sigma
#Time Window
tspan = (0.0,T)
#Definition of the SDE
prob = SDEProblem(f1,g1_,u0,tspan)
rate(u,p,t) = (lam*T)
affect!(integrator) = (integrator.u = integrator.u+((rand()<p1) ? quantile_exp(lamp,rand()) : -quantile_exp(lamm,rand())))
jump = ConstantRateJump(rate,affect!)
jump_prob = JumpProblem(prob,Direct(),jump)
Model = MonteCarloProblem(jump_prob)


FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)

@show FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData);

@test abs(FwdPrice-99.41332633109904)<toll
@test abs(EuPrice-10.347332240535199)<toll
@test abs(BarrierPrice-8.860123655599818)<toll
@test abs(AsianPrice1-5.81798437145069)<toll