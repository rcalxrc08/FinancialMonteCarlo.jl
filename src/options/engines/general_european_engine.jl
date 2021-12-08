
function payoff(S::AbstractMatrix{num}, euPayoff::EuropeanPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(euPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = euPayoff.T
    NStep = mcBaseData.Nstep
    index1 = round(UInt, T / T1 * NStep) + 1
    @views ST = S[:, index1]
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, euPayoff)
    payoff2 = phi.(ST)

    return payoff2 * exp(-integral(r, T))
end
