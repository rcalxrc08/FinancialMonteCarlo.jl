using DelimitedFiles

function simulate(mcProcess::BrownianMotion,spotData::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: PrescribedMC}
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
	Z=DelimitedFiles.readdlm("C:\\Users\\Nicola\\.julia\\dev\\FinancialMonteCarlo\\src\\models\\matrix.txt")[1:Nsim,1:Nstep];
	view(X,:,1).=isDualZero;
	@inbounds for i=1:Nsim
		@inbounds for j=1:Nstep
			x_i_j=@views X[i,j];
			z_i_j=@views Z[i,j];
			X[i,j+1]=x_i_j+mean_bm+stddev_bm*z_i_j;
		end
	end

	return X;

end
