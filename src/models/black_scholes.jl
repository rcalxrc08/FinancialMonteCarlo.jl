type BlackScholesProcess{num<:Number}<:ItoProcess
	sigma::num
	function BlackScholesProcess(sigma::num) where {num <: Number}
        if sigma <= 0.0
            error("Volatility must be positive")
        else
            return new{num}(sigma)
        end
    end
end

export BlackScholesProcess;

function simulate(mcProcess::BlackScholesProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,monteCarloMode::MonteCarloMode=standard)
	if T<=0.0
		error("Final time must be positive");
	end
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	sigma_gbm=mcProcess.sigma;
	mu_gbm=r-d;
	GeomData=MonteCarloBaseData(mcBaseData.Nsim,mcBaseData.Nstep)
	
	S=S0.*simulate(GeometricBrownianMotion(sigma_gbm,mu_gbm),spotData,GeomData,T,monteCarloMode)
	
	return S;
	
end