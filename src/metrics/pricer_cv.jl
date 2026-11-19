# The expectation of the discounted control: the price given to ControlVariates, the closed form of a forward, or the
# mean of the pilot, which cancels most of the variance reduction (#75)
control_price(price::Number, variate, variate_values, mcProcess::BaseProcess, T) = price
control_price(::Nothing, ::AbstractPayoff, variate_values, ::BaseProcess, T) = mean(variate_values)
control_price(::Nothing, ::Forward, variate_values, mcProcess::BaseProcess, T) = mcProcess.underlying.S0 * exp(-integral(mcProcess.underlying.d, T))

function pricer_control_variates(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff)
    set_seed!(mcConfig)
    variate_handl = mcConfig.sampling
    variate_conf = variate_handl.conf_variate
    variate_payoff = variate_handl.variate
    set_seed!(variate_conf)
    ChainRulesCore.@ignore_derivatives check_argument(maturity(abstractPayoff) == variate_payoff.T, "The control variate must have the maturity of the payoff")
    T = maturity(abstractPayoff)
    S_var = simulate(mcProcess, rfCurve, variate_conf, T)
    variate_values = payoff(S_var, variate_payoff, rfCurve, variate_conf)
    option_values = payoff(S_var, abstractPayoff, rfCurve, variate_conf)
    # the optimal coefficient, estimated on the pilot
    c = -cov(variate_values, option_values) / var(variate_values)
    price_var = control_price(variate_handl.price, variate_payoff, variate_values, mcProcess, T)
    mcConfig_mod = @set mcConfig.sampling = variate_handl.sampling
    prices = pricer(mcProcess, rfCurve, mcConfig_mod, [abstractPayoff, variate_payoff])
    @muladd price = prices[1] + c * (prices[2] - price_var)

    return price
end

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:ControlVariates}, abstractPayoff::AbstractPayoff)
    return pricer_control_variates(mcProcess, rfCurve, mcConfig, abstractPayoff)
end
