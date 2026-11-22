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

# The price of the spot is S0 whatever the configuration: a method for each configuration type with a pricer of its own,
# so that none is ambiguous (#41)
for config_type in (MonteCarloConfiguration, BatchedConfig, MonteCarloConfiguration{<:ControlVariates}, MonteCarloConfiguration{<:ControlVariates, <:AbstractBatchedMode})
    @eval pricer(mcProcess::BaseProcess, ::AbstractZeroRateCurve, ::$config_type, ::Spot) = mcProcess.underlying.S0
end

# No maturity: an Int8 zero, which does not promote the maturities of the other payoffs of a position (#39)
maturity(::Spot) = Int8(0)
