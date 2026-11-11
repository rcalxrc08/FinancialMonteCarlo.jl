
function payoff(S::AbstractMatrix{num}, payoff_::PathDependentPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = payoff_.T
    Nsim = mcBaseData.Nsim
    index1 = maturity_index(T, T1, mcBaseData)

    @inbounds f(S::abstractArray) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number} = payout(S, payoff_)
    @inbounds payoff2 = [f(view(S, i, 1:index1)) for i = 1:Nsim]

    return payoff2 * exp(-integral(r, T))
end

# The same payoffs from reductions of the rows of S, with no path indexed one by one: the minimum and maximum of each path
# for the barrier options, its mean for the Asian ones. The CUDA extension prices them so on the device; other path
# dependent payoffs collect the paths to the host.
function payoff_by_rows(S::AbstractMatrix, payoff_::BarrierPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number = maturity(payoff_))
    T = payoff_.T
    index1 = maturity_index(T, T1, mcBaseData)
    S_ = view(S, :, 1:index1)
    lo = vec(minimum(S_; dims = 2))
    hi = vec(maximum(S_; dims = 2))
    return barrier_payout.(values_at_maturity(S, T, T1, mcBaseData), lo, hi, payoff_) .* exp(-integral(rfCurve.r, T))
end
function payoff_by_rows(S::AbstractMatrix, payoff_::AsianPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number = maturity(payoff_))
    T = payoff_.T
    index1 = maturity_index(T, T1, mcBaseData)
    A = vec(mean(view(S, :, 1:index1); dims = 2))
    return asian_payout.(values_at_maturity(S, T, T1, mcBaseData), A, payoff_) .* exp(-integral(rfCurve.r, T))
end
