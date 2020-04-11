using Test, FinancialMonteCarlo,DualNumbers;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

underlying_name="ENJ"

Nsim=10000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8

rfCurve=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

display(Model)

port_=EUData*2.0+FwdData-FwdData/2.0-FwdData;
port_=port_*2.0;
port_=2.0*port_;
display(port_)
print(port_)
display(FwdData)
display(-FwdData)
print(FwdData)
print(-FwdData)
display(Forward(dual(T,1.0)))
display(-Forward(dual(T,1.0)))
print(Forward(dual(T,1.0)))
print(-Forward(dual(T,1.0)))

display(Forward(dual(T,1.0))+EUData)
display(-Forward(dual(T,1.0))+EUData)
print(Forward(dual(T,1.0))+EUData)
print(-Forward(dual(T,1.0))+EUData)

display(Forward(dual(T,1.0))-EUData)
display(-Forward(dual(T,1.0))-EUData)
print(Forward(dual(T,1.0))-EUData)
print(-Forward(dual(T,1.0))-EUData)
fwd=Forward(dual(T,1.0))*1.0
fwd=1.0*Forward(dual(T,1.0))
ok=FwdData-fwd;
display(fwd)
print(fwd)
fwd=1.0*fwd
fwd=fwd*1.0
fwd=fwd/1.0
fwd=-fwd
fwd=Forward(dual(T,1.0))*dual(1.0,1.0)+EUData*2.0+AMData
display(fwd)
print(fwd)
fwd=-Forward(dual(T,1.0))*dual(1.0,1.0)+EUData*2.0+AMData
display(fwd)
print(fwd)

mktdataset=underlying_name→Model
portfolio_=[FwdData;EUData;AMData;BarrierData;AsianFloatingStrikeData;AsianFixedStrikeData;Spot()];
portfolio=underlying_name→FwdData
portfolio+=underlying_name→EUData
portfolio+=underlying_name→AMData
portfolio+=underlying_name→Spot()
portfolio+=underlying_name→BarrierData
portfolio+=underlying_name→AsianFloatingStrikeData
portfolio+=underlying_name→AsianFixedStrikeData

display(portfolio)


price_mkt=pricer(mktdataset,rfCurve,mc,portfolio)
price_old= sum(pricer(Model,rfCurve,mc,portfolio_))


@test abs(price_mkt-price_old)<1e-8