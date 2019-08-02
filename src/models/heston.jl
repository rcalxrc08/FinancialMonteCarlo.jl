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
mutable struct HestonProcess{num <: Number, num1 <: Number , num2 <: Number , num3 <: Number , num4 <: Number, num5 <:Number}<:ItoProcess
	σ::num
	σ_zero::num1
	λ::num2
	κ::num3
	ρ::num4
	θ::num5
	function HestonProcess(σ::num,σ_zero::num1,λ::num2,	κ::num3,ρ::num4,θ::num5) where {num <: Number, num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number,num5 <: Number}
        if σ_zero<=0.0
			error("initial volatility must be positive");
		elseif σ<=0.0
			error("volatility of volatility must be positive");
		elseif abs(κ+λ)<=1e-14
			error("unfeasible parameters");
		elseif !(-1.0<=ρ<=1.0)
			error("ρ must be a correlation");
        else
            return new{num,num1,num2,num3,num4,num5}(σ,σ_zero,λ,κ,ρ,θ)
        end
    end
end

export HestonProcess;

function simulate(mcProcess::HestonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3},T::numb,parallelMode::SerialMode=SerialMode()) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC}
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	σ_zero=mcProcess.σ_zero;
	λ1=mcProcess.λ;
	κ=mcProcess.κ;
	ρ=mcProcess.ρ;
	θ=mcProcess.θ;
	if T<=0.0
		error("Final time must be positive");
	end

	####Simulation
	## Simulate
	κ_s=κ+λ1;
	θ_s=κ*θ/(κ+λ1);

	dt=T/Nstep
	isDualZero=S0*T*r*σ_zero*κ*θ*λ1*σ*ρ*0.0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	v_m=[σ_zero+isDualZero for _ in 1:Nsim];
	for j in 1:Nstep
		e1=randn(Nsim);
		e2=randn(Nsim);
		e2=e1.*ρ.+e2.*sqrt(1-ρ*ρ);
		X[:,j+1]=X[:,j]+((r-d).-0.5.*max.(v_m,0)).*dt+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
		v_m+=κ_s.*(θ_s.-max.(v_m,0)).*dt+σ.*sqrt.(max.(v_m,0)).*sqrt(dt).*e2;
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end

function simulate(mcProcess::HestonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3},T::numb,parallelMode::SerialMode=SerialMode()) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC}
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	σ_zero=mcProcess.σ_zero;
	λ1=mcProcess.λ;
	κ=mcProcess.κ;
	ρ=mcProcess.ρ;
	θ=mcProcess.θ;
	if T<=0.0
		error("Final time must be positive");
	end

	####Simulation
	## Simulate
	κ_s=κ+λ1;
	θ_s=κ*θ/(κ+λ1);

	dt=T/Nstep
	isDualZero=S0*T*r*σ_zero*κ*θ*λ1*σ*ρ*0.0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	v_m=[σ_zero+isDualZero for _ in 1:Nsim];
	Nsim_2=Int(floor(Nsim/2))
	if Nsim_2*2!=Nsim
		error("Antithetic support only odd number of simulations")
	end
	for j in 1:Nstep
		e1=randn(Nsim_2);
		e2=randn(Nsim_2);
		e1=[e1;-e1];
		e2=[e2;-e2];
		e2=e1.*ρ.+e2.*sqrt(1-ρ*ρ);
		X[:,j+1]=X[:,j]+((r-d).-0.5.*max.(v_m,0)).*dt+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
		v_m+=κ_s.*(θ_s.-max.(v_m,0)).*dt+(σ*sqrt(dt)).*sqrt.(max.(v_m,0)).*e2;
	end
	
	## Conclude
	S=S0.*exp.(X);
	return S;

end