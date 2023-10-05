"""
Struct for Black Scholes Process

		bsProcess=BlackScholesProcess(σ::num1) where {num1 <: Number}
	
Where:

		σ	=	volatility of the process.
"""
struct BlackScholesProcess{num <: Number, abstrUnderlying <: AbstractUnderlying} <: ItoProcess{num}
    σ::num
    underlying::abstrUnderlying
    function BlackScholesProcess(σ::num, underlying::abstrUnderlying) where {num <: Number, abstrUnderlying <: AbstractUnderlying}
        @assert σ > 0 "Volatility must be positive"
        return new{num, abstrUnderlying}(σ, underlying)
    end
end

export BlackScholesProcess;

function simulate!(S, mcProcess::BlackScholesProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    @assert T > 0.0
    r = rfCurve.r
    d = dividend(mcProcess)
    σ = mcProcess.σ
    μ = r - d

    simulate!(S, GeometricBrownianMotion(σ, μ, mcProcess.underlying.S0), mcBaseData, T)

    nothing
end

# function simulate_path!(S, mcProcess::BlackScholesProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
#     @assert T > 0.0
#     r = rfCurve.r
#     d = dividend(mcProcess)
#     σ = mcProcess.σ
#     μ = r - d

#     simulate_path!(S, GeometricBrownianMotion(σ, μ, mcProcess.underlying.S0), mcBaseData, T)

#     nothing
# end
