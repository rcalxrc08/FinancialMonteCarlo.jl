
function simulate(mcProcess::HestonProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode,parallelMode::CudaMode) where {numb<:Number}
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
	if monteCarloMode==antithetic
		Nsim_2=Int(floor(Nsim/2))
		if Nsim_2*2!=Nsim
			error("Antithetic support only odd number of simulations")
		end
		ρ_f=Float32(ρ);
		κ_s_f=Float32(κ_s);
		θ_s_f=Float32(θ_s);
		r_d_f=Float32(r-d)
		dt_f=Float32(dt)
		σ_f=Float32(σ)
		X=CuMatrix{typeof(Float32(isDualZero+σ_zero))}(undef,Nsim,Nstep+1);
		X[:,1].=Float32(isDualZero);
		v_m=CuArray{typeof(σ_zero+isDualZero)}(undef,Nsim)
		for j in 1:Nstep
			e1=randn(Float32,Nsim_2);
			e2=randn(Float32,Nsim_2);
			e1=cu([e1;-e1]);
			e2=cu([e2;-e2]);
			e2=e1.*ρ_f.+e2.*sqrt(1-ρ_f*ρ_f);
			X[:,j+1]=X[:,j]+((r_d_f).-0.5.*max.(v_m,0)).*dt_f+sqrt.(max.(v_m,0)).*sqrt(dt_f).*e1;
			v_m+=κ_s_f.*(θ_s_f.-max.(v_m,0)).*dt_f+σ_f.*sqrt.(max.(v_m,0)).*sqrt(dt_f).*e2;
		end
		## Conclude
		S=S0.*exp.(Matrix(X));
		return S;
	else
		ρ_f=Float32(ρ);
		κ_s_f=Float32(κ_s);
		θ_s_f=Float32(θ_s);
		r_d_f=Float32(r-d)
		dt_f=Float32(dt)
		σ_f=Float32(σ)
		X=CuMatrix{typeof(Float32(isDualZero+σ_zero))}(undef,Nsim,Nstep+1);
		X[:,1].=Float32(isDualZero);
		v_m=CuArray{typeof(σ_zero+isDualZero)}(undef,Nsim)
		for j in 1:Nstep
			e1=cu(randn(Float32,Nsim));
			e2=cu(randn(Float32,Nsim));
			e2=e1.*ρ_f.+e2.*sqrt(1-ρ_f*ρ_f);
			X[:,j+1]=X[:,j]+((r_d_f).-0.5.*max.(v_m,0)).*dt_f+sqrt.(max.(v_m,0)).*sqrt(dt_f).*e1;
			v_m+=κ_s_f.*(θ_s_f.-max.(v_m,0)).*dt_f+σ_f.*sqrt.(max.(v_m,0)).*sqrt(dt_f).*e2;
		end
		## Conclude
		S=S0.*exp.(Matrix(X));
		return S;
	end

end
