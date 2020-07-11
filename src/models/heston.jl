"""
Struct for Heston Process

		hstProcess=HestonProcess(σ::num1,σ_zero::num2,λ::num3,κ::num4,ρ::num5,θ::num6) where {num1,num2,num3,num4,num5,num6 <: Number}
	
Where:\n
		σ	=	volatility of volatility of the process.
		σ_zero	= initial volatility of the process.
		λ	=	??volatility of the process.
		κ	=	??volatility of the process.
		ρ	=	??volatility of the process.
		θ	=	??drift of the process.
"""
mutable struct HestonProcess{num <: Number, num1 <: Number , num2 <: Number , num3 <: Number , num4 <: Number, num5 <:Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ItoProcess{numtype}
	σ::num
	σ_zero::num1
	λ::num2
	κ::num3
	ρ::num4
	θ::num5
	underlying::abstrUnderlying
	function HestonProcess(σ::num,σ_zero::num1,λ::num2,	κ::num3,ρ::num4,θ::num5,underlying::abstrUnderlying) where {num <: Number, num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number,num5 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ_zero<=0.0
			error("initial volatility must be positive");
		elseif σ<=0.0
			error("volatility of volatility must be positive");
		elseif abs(κ+λ)<=1e-14
			error("unfeasible parameters");
		elseif !(-1.0<=ρ<=1.0)
			error("ρ must be a correlation");
        else
			zero_typed=zero(num)+zero(num1)+zero(num2)+zero(num3)+zero(num4)+zero(num5);
            return new{num,num1,num2,num3,num4,num5,abstrUnderlying,typeof(zero_typed)}(σ,σ_zero,λ,κ,ρ,θ,underlying)
        end
    end
end

export HestonProcess;

function simulate!(X,mcProcess::HestonProcess,rfCurve::ZeroRate,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: StandardMC, type5 <: Random.AbstractRNG}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	σ_zero=mcProcess.σ_zero;
	λ1=mcProcess.λ;
	κ=mcProcess.κ;
	ρ=mcProcess.ρ;
	θ=mcProcess.θ;
	@assert T>0.0

	####Simulation
	## Simulate
	κ_s=κ+λ1;
	θ_s=κ*θ/(κ+λ1);

	dt=T/Nstep
	isDualZero=T*r*σ_zero*θ_s*κ_s*σ*ρ*0.0*S0;
	#X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	view(X,:,1).=isDualZero;
	for i in 1:Nsim
		v_m=σ_zero^2;
		for j in 1:Nstep
			e1=randn(mcBaseData.rng);
			e2=e1*ρ+randn(mcBaseData.rng)*sqrt(1-ρ*ρ);
			@views X[i,j+1]=X[i,j]+((r-d)-0.5*v_m)*dt+sqrt(v_m)*sqrt(dt)*e1;
			v_m+=κ_s*(θ_s-v_m)*dt+σ*sqrt(v_m)*sqrt(dt)*e2;
			#when v_m = 0.0, the derivative becomes NaN
			v_m=max(v_m,isDualZero);
		end
	end
	## Conclude
	X.=S0.*exp.(X);
	return;

end

function simulate!(X,mcProcess::HestonProcess,rfCurve::ZeroRate,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	σ_zero=mcProcess.σ_zero;
	λ1=mcProcess.λ;
	κ=mcProcess.κ;
	ρ=mcProcess.ρ;
	θ=mcProcess.θ;
	@assert T>0.0

	####Simulation
	## Simulate
	κ_s=κ+λ1;
	θ_s=κ*θ/(κ+λ1);

	dt=T/Nstep
	isDualZero=T*r*σ_zero*κ*θ*λ1*σ*ρ*0.0;
	#X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	view(X,:,1).=isDualZero;
	for i in 1:div(Nsim,2)
		v_m_1=σ_zero^2;
		v_m_2=σ_zero^2;
		for j in 1:Nstep
			e1=randn(mcBaseData.rng);
			e2= -(e1*ρ+randn(mcBaseData.rng)*sqrt(1-ρ*ρ));
			@views X[2*i-1,j+1]=X[2*i-1,j]+((r-d)-0.5*v_m_1)*dt+sqrt(v_m_1)*sqrt(dt)*e1;
			@views X[2*i,j+1]=X[2*i,j]+((r-d)-0.5*v_m_2)*dt+sqrt(v_m_2)*sqrt(dt)*(-e1);
			v_m_1+=κ_s*(θ_s-v_m_1)*dt+σ*sqrt(v_m_1)*sqrt(dt)*e2;
			v_m_2+=κ_s*(θ_s-v_m_2)*dt+σ*sqrt(v_m_2)*sqrt(dt)*(-e2);
			v_m_1=max(v_m_1,isDualZero);
			v_m_2=max(v_m_2,isDualZero);
		end
	end
	## Conclude
	X.=S0.*exp.(X);
	return;

end