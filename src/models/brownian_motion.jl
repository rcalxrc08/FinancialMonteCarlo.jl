"""
Struct for Brownian Motion

		bmProcess=BrownianMotion(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct BrownianMotion{num <: Number, num1 <: Number, num2 <: Number , num3 <: Number} <: ItoProcess
	σ::num
	μ::num1
	underlying::Underlying{num2,num3}
	function BrownianMotion(σ::num,μ::num1,underlying::Underlying{num2,num3}) where {num <: Number, num1 <: Number, num2 <: Number , num3 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,num2,num3}(σ,μ,underlying)
        end
    end
end

export BrownianMotion;

function simulate(mcProcess::BrownianMotion,spotData::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: StandardMC}
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
	isDualZero=mean_bm*stddev_bm*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	view(X,:,1).=isDualZero;
	@inbounds for i=1:Nsim
		@inbounds for j=1:Nstep
			x_i_j=@views X[i,j];
			X[i,j+1]=x_i_j+mean_bm+stddev_bm*randn(mcBaseData.rng);
		end
	end

	return X;

end


function simulate(mcProcess::BrownianMotion,spotData::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AntitheticMC}
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
	isDualZero=mean_bm*stddev_bm*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	Nsim_2=div(Nsim,2)

	@inbounds for j in 1:Nstep
		@inbounds for i in 1:Nsim_2
			Z=stddev_bm.*randn(mcBaseData.rng);
			X[2*i-1,j+1]=X[2*i-1,j].+mean_bm.+Z;
			X[2*i,j+1]  =X[2*i,j]  .+mean_bm.-Z;
		end
	end

	return X;

end