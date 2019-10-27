function pricer_macro_2(model_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: ControlVariates ,  <: BaseMode},abstractPayoff::AbstractPayoff)

			set_seed(mcConfig)
			variate_handl=mcConfig.monteCarloMethod
			variate_conf=variate_handl.conf_variate;
			variate_payoff=variate_handl.variate;
			@assert abstractPayoff.T==variate_payoff.T
			T=abstractPayoff.T;
			S_var=simulate(mcProcess,spotData,variate_conf,T)
			Payoff_var=payoff(S_var,variate_payoff,spotData);
			Payoff_opt_var=payoff(S_var,abstractPayoff,spotData);
			c=-cov(Payoff_var,Payoff_opt_var)/var(Payoff_var);
			price_var=mean(Payoff_var);
			mcConfig_mod=MonteCarloConfiguration(mcConfig.Nsim,mcConfig.Nstep,variate_conf.monteCarloMethod,mcConfig.parallelMode,mcConfig.seed*2)
			#END OF VARIATE SECTION
			Prices=pricer(mcProcess,spotData,mcConfig_mod,[abstractPayoff,variate_payoff]);
			Price=Prices[1]+c*(Prices[2]-price_var);
			return Price;
		end
	end
end
pricer_macro_2(BaseProcess)


function pricer_macro_3(model_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: ControlVariates{Forward{num}} ,  <: BaseMode},abstractPayoff::AbstractPayoff) where { num <: Number }
			@show "special"
			set_seed(mcConfig)
			variate_handl=mcConfig.monteCarloMethod
			variate_conf=variate_handl.conf_variate;
			variate_payoff=variate_handl.variate;
			@assert abstractPayoff.T==variate_payoff.T
			T=abstractPayoff.T;
			eps_1=0.01
			r=spotData.r;
			S_var=simulate(mcProcess,spotData,variate_conf,T)
			Payoff_var=payoff(S_var,variate_payoff,spotData);
			Payoff_opt_var=payoff(S_var,abstractPayoff,spotData);
			c=cov(Payoff_var,Payoff_opt_var)/var(Payoff_var);
			@show mean(Payoff_var)
			@show c
			price_var=spotData.S0;
			mcConfig_mod=MonteCarloConfiguration(mcConfig.Nsim,mcConfig.Nstep,variate_conf.monteCarloMethod,mcConfig.parallelMode,mcConfig.seed+3)
			#END OF VARIATE SECTION
			Prices=pricer(mcProcess,spotData,mcConfig_mod,[abstractPayoff,variate_payoff]);
			@show Prices[2]
			Price=Prices[1]-eps_1*c*(Prices[2]-price_var);
			@show Prices[1]+c*(Prices[2]-price_var);
			return Price;
		end
	end
end
pricer_macro_3(BaseProcess)