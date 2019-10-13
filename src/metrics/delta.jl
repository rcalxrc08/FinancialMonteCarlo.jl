
function delta_macro(model_type)
	@eval begin
		"""
		General Interface for Pricing

				Delta=delta(mcProcess,spotData,mcBaseData,payoff_,parallelMode=SerialMode())
			
		Where:\n
				mcProcess          = Process to be simulated.
				spotData  = Datas of the Spot.
				mcBaseData = Basic properties of MonteCarlo simulation
				payoff_ = Payoff(s) to be priced
				

				Delta     = Delta of the derivative

		"""	
		function delta(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,dS0::Real=1e-7)

			Random.seed!(0)
			Price=pricer(mcProcess,spotData,mcConfig,abstractPayoff);
			Random.seed!(0)
			spotData_1=equitySpotData(spotData.S0+dS0,spotData.r,spotData.d);
			PriceUp=pricer(mcProcess,spotData_1,mcConfig,abstractPayoff);
			Delta=(PriceUp-Price)/dS0;

			return Delta;
		end
	end
end

delta_macro(BaseProcess)

function delta_macro_array(model_type)
	@eval function delta(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},dS0::Real=1e-7)
			Random.seed!(0)
			Prices=pricer(mcProcess,spotData,mcConfig,abstractPayoffs);
			Random.seed!(0)
			spotData_1=equitySpotData(spotData.S0+dS0,spotData.r,spotData.d);
			PricesUp=pricer(mcProcess,spotData_1,mcConfig,abstractPayoffs);
			Delta=(PricesUp.-Prices)./dS0;
		
		return Delta;
	end
end

delta_macro_array(BaseProcess)