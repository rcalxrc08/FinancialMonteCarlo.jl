
type GeometricBrownianMotion{num,num1<:Number}<:ItoProcess
	σ::num
	μ::num1
	function GeometricBrownianMotion(σ::num,μ::num1) where {num,num1 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1}(σ,μ)
        end
    end
end

export GeometricBrownianMotion;

function simulate(mcProcess::GeometricBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	if T<=0.0
		error("Final time must be positive");
	end
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-σ_gbm^2/2;
	const dictGBM=Dict{String,Number}("σ"=>σ_gbm, "μ" => μ_bm)
	BrownianData=MonteCarloConfiguration(mcBaseData.Nsim,mcBaseData.Nstep)
	X=simulate(BrownianMotion(σ_gbm,μ_bm),spotData,BrownianData,T,monteCarloMode)
	S=exp.(X);
	return S;
end