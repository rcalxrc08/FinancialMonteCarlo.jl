function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:ControlVariates, <:BaseMode}, abstractPayoff::AbstractPayoff)
    set_seed(mcConfig)
    variate_handl = mcConfig.monteCarloMethod
    variate_conf = variate_handl.conf_variate
    variate_payoff = variate_handl.variate
    @assert maturity(abstractPayoff) == variate_payoff.T
    T = maturity(abstractPayoff)
    S_var = simulate(mcProcess, rfCurve, variate_conf, T)
    Payoff_var = payoff(S_var, variate_payoff, rfCurve)
    Payoff_opt_var = payoff(S_var, abstractPayoff, rfCurve)
    #c=-cov(Payoff_var,Payoff_opt_var)/var(Payoff_var)[1];
    c = -collect(cov(Payoff_var, Payoff_opt_var) / var(Payoff_var))[1]
    price_var = mean(Payoff_var)
    mcConfig_mod = MonteCarloConfiguration(mcConfig.Nsim, mcConfig.Nstep, variate_conf.monteCarloMethod, mcConfig.parallelMode, mcConfig.seed * 2)
    #END OF VARIATE SECTION
    Prices = pricer(mcProcess, rfCurve, mcConfig_mod, [abstractPayoff, variate_payoff])
    Price = Prices[1] + c * (Prices[2] - price_var)

    return Price
end

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:ControlVariates{Forward{num}}, <:BaseMode}, abstractPayoff::AbstractPayoff) where {num <: Number}
    set_seed(mcConfig)
    variate_handl = mcConfig.monteCarloMethod
    variate_conf = variate_handl.conf_variate
    variate_payoff = variate_handl.variate
    @assert maturity(abstractPayoff) == variate_payoff.T
    T = maturity(abstractPayoff)
    S_var = simulate(mcProcess, rfCurve, variate_conf, T)
    Payoff_var = payoff(S_var, variate_payoff, rfCurve)
    Payoff_opt_var = payoff(S_var, abstractPayoff, rfCurve)
    c = -collect(cov(Payoff_var, Payoff_opt_var) / var(Payoff_var))[1]
    price_var = mcProcess.underlying.S0 * exp(-integral(mcProcess.underlying.d, T))
    mcConfig_mod = MonteCarloConfiguration(mcConfig.Nsim, mcConfig.Nstep, variate_conf.monteCarloMethod, mcConfig.parallelMode, mcConfig.seed + 3)
    #END OF VARIATE SECTION
    Prices = pricer(mcProcess, rfCurve, mcConfig_mod, [abstractPayoff, variate_payoff])
    Price = Prices[1] - c * (Prices[2] - price_var)
    return Price
end
