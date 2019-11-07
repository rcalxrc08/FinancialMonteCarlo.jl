using FinancialMonteCarlo
@show "Test Structs Building"
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
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8;

spotData1=ZeroRateCurve(r);
@test_throws(ErrorException,MonteCarloConfiguration(-Nsim,Nstep))
@test_throws(ErrorException,MonteCarloConfiguration(Nsim,-Nstep))
#####################################################################
#Payoff Structs Test: Negative TTM
@test_throws(ErrorException,Forward(Tneg));
@test_throws(ErrorException,EuropeanOption(Tneg,K));
@test_throws(ErrorException,AmericanOption(Tneg,K));
@test_throws(ErrorException,BinaryAmericanOption(Tneg,K));
@test_throws(ErrorException,BinaryEuropeanOption(Tneg,K));
@test_throws(ErrorException,BarrierOptionDownOut(Tneg,K,D));
@test_throws(ErrorException,BarrierOptionUpOut(Tneg,K,D));
@test_throws(ErrorException,BarrierOptionUpIn(Tneg,K,D));
@test_throws(ErrorException,BarrierOptionDownIn(Tneg,K,D));
@test_throws(ErrorException,AsianFloatingStrikeOption(Tneg));
@test_throws(ErrorException,AsianFixedStrikeOption(Tneg,K));
@test_throws(ErrorException,DoubleBarrierOption(Tneg,K,K*10,D));

# Negative Strike and Barriers
@test_throws(ErrorException,EuropeanOption(T,Kneg));
@test_throws(ErrorException,AmericanOption(T,Kneg));
@test_throws(ErrorException,BinaryEuropeanOption(T,Kneg));
@test_throws(ErrorException,BinaryAmericanOption(T,Kneg));
@test_throws(ErrorException,AsianFixedStrikeOption(T,Kneg));
@test_throws(ErrorException,BarrierOptionDownOut(T,Kneg,D));
@test_throws(ErrorException,BarrierOptionDownOut(T,K,Kneg));
@test_throws(ErrorException,BarrierOptionDownIn(T,Kneg,D));
@test_throws(ErrorException,BarrierOptionDownIn(T,K,Kneg));
@test_throws(ErrorException,BarrierOptionUpIn(T,Kneg,D));
@test_throws(ErrorException,BarrierOptionUpIn(T,K,Kneg));
@test_throws(ErrorException,BarrierOptionUpOut(T,Kneg,D));
@test_throws(ErrorException,BarrierOptionUpOut(T,K,Kneg));
@test_throws(ErrorException,DoubleBarrierOption(T,K,K,Kneg));
@test_throws(ErrorException,DoubleBarrierOption(T,Kneg,K,K));
@test_throws(ErrorException,DoubleBarrierOption(T,K,Kneg,K));