control_price(::AbstractPayoff, Payoff_var, ::BaseProcess, T) = mean(Payoff_var)
control_price(::Forward, Payoff_var, mcProcess::BaseProcess, T) = mcProcess.underlying.S0 * exp(-integral(mcProcess.underlying.d, T))

function pricer_control_variates(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff)
    set_seed!(mcConfig)
    variate_handl = mcConfig.sampling
    variate_conf = variate_handl.conf_variate
    variate_payoff = variate_handl.variate
    set_seed!(variate_conf)
    ChainRulesCore.@ignore_derivatives @assert maturity(abstractPayoff) == variate_payoff.T
    T = maturity(abstractPayoff)
    S_var = simulate(mcProcess, rfCurve, variate_conf, T)
    Payoff_var = payoff(S_var, variate_payoff, rfCurve, variate_conf)
    Payoff_opt_var = payoff(S_var, abstractPayoff, rfCurve, variate_conf)
    #c=-cov(Payoff_var,Payoff_opt_var)/var(Payoff_var)[1];
    c = -collect(cov(Payoff_var, Payoff_opt_var) / var(Payoff_var))[1]
    price_var = control_price(variate_payoff, Payoff_var, mcProcess, T)
    mcConfig_mod = @set mcConfig.sampling = variate_handl.sampling
    #END OF VARIATE SECTION
    Prices = pricer(mcProcess, rfCurve, mcConfig_mod, [abstractPayoff, variate_payoff])
    Price = Prices[1] + c * (Prices[2] - price_var)

    return Price
end

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:ControlVariates}, abstractPayoff::AbstractPayoff)
    return pricer_control_variates(mcProcess, rfCurve, mcConfig, abstractPayoff)
end
