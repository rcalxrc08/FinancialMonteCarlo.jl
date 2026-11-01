
values_at_maturity(S::AbstractMatrix, T, T1, mcBaseData::AbstractMonteCarloConfiguration) = @views S[:, maturity_index(T, T1, mcBaseData)]

function payoff(S::AbstractMatrix{num}, euPayoff::EuropeanPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(euPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = euPayoff.T
    ST = values_at_maturity(S, T, T1, mcBaseData)
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, euPayoff)
    payoff2 = @. phi(ST)

    return payoff2 * exp(-integral(r, T))
end
