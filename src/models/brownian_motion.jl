
mutable struct BrownianMotion{num,num1<:Number}<:ItoProcess
	σ::num
	μ::num1
	function BrownianMotion(σ::num,μ::num1) where {num,num1 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1}(σ,μ)
        end
    end
end

export BrownianMotion;

Float32(p::Dual{Float64})=dual(Float32(p.value),Float32(p.epsilon))

function simulate(mcProcess::BrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard) where {numb<:Number}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;

	σ=mcProcess.σ;
	μ=mcProcess.μ;
	if T<=0.0
		error("Final time must be positive");
	end
	dt=T/Nstep
	mean_bm=μ*dt
	stddev_bm=σ*sqrt(dt)
	if monteCarloMode==antithetic
		isDualZero=mean_bm*stddev_bm*0.0;
		X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
		X[:,1].=isDualZero;
		Nsim_2=Int(floor(Nsim/2))
		if Nsim_2*2!=Nsim
			error("Antithetic support only odd number of simulations")
		end
		for j in 1:Nstep
			Z=randn(Nsim_2);
			Z=[Z;-Z];
			X[:,j+1]=X[:,j].+mean_bm.+stddev_bm.*Z;
		end
		return X;
	elseif monteCarloMode==parallel_cuda_gpu
		mean_bm_f=Float32(mean_bm);
		stddev_bm_f=Float32(stddev_bm);
		isDualZero=mean_bm_f*stddev_bm_f*zero(Float32);
		X_cu=CuMatrix{typeof(isDualZero)}(Nsim,Nstep+1);
		for i=1:Nstep
			X_cu[:,i+1]=X_cu[:,i]+(mean_bm_f.+stddev_bm_f.*cu(randn(Float32,Nsim)));
		end
		X=Matrix(X_cu);
		return X;
	else
		isDualZero=mean_bm*stddev_bm*0.0;
		X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
		X[:,1].=isDualZero;
		for i=1:Nsim
			for j=1:Nstep
				X[i,j+1]=X[i,j].+mean_bm.+stddev_bm.*randn();
			end
		end
		return X;
	end

end
