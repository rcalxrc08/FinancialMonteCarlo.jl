
mutable struct HestonProcess{num,num1,num2,num3,num4,num5<:Number}<:ItoProcess
	σ::num
	σ_zero::num1
	λ::num2
	κ::num3
	ρ::num4
	θ::num5
	function HestonProcess(σ::num,σ_zero::num1,λ::num2,	κ::num3,ρ::num4,θ::num5) where {num,num1,num2,num3,num4,num5 <: Number}
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

function simulate(mcProcess::HestonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard) where {numb<:Number}
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
	if monteCarloMode!=antithetic
		v_m=[σ_zero+isDualZero for _ in 1:Nsim];
		for j in 1:Nstep
			e1=randn(Nsim);
			e2=randn(Nsim);
			e2=e1.*ρ.+e2.*sqrt(1-ρ*ρ);
			X[:,j+1]=X[:,j]+((r-d).-0.5.*max.(v_m,0)).*dt+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
			v_m+=κ_s.*(θ_s.-max.(v_m,0)).*dt+σ.*sqrt.(max.(v_m,0)).*sqrt(dt).*e2;
		end
	else
		v_m=[σ_zero+isDualZero for _ in 1:Nsim];
		NsimAnti=Int(floor(Nsim/2))
		for j in 1:Nstep
			e1=randn(NsimAnti);
			e2=randn(NsimAnti);
			e1=[e1;-e1];
			e2=[e2;-e2];
			e2=e1.*ρ.+e2.*sqrt(1-ρ*ρ);
			X[:,j+1]=X[:,j]+((r-d).-0.5.*max.(v_m,0)).*dt+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
			v_m+=κ_s.*(θ_s.-max.(v_m,0)).*dt+(σ*sqrt(dt)).*sqrt.(max.(v_m,0)).*e2;
		end
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end
