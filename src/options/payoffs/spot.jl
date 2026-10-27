"""
    Spot()

The spot price of the underlying today: its price is `S0`, whatever the model.
Mostly useful inside positions and portfolios.
"""
struct Spot <: EuropeanPayoff{Float16}
    function Spot()
        return new()
    end
end

export Spot;

function payoff(S::AbstractMatrix{num}, optionData::Spot, ::ZeroRate, ::AbstractMonteCarloConfiguration, T1::num2 = maturity(optionData)) where {num <: Number, num2 <: Number}
    S0_vec = S[:, 1]

    return S0_vec
end

function pricer(mcProcess::BaseProcess, ::AbstractZeroRateCurve, ::MonteCarloConfiguration, ::Spot)
    return mcProcess.underlying.S0
end

function maturity(::Spot)
    return 0.0
end
