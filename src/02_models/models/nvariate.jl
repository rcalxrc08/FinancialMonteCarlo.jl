"""
    GaussianCopulaNVariateProcess(ρ, models...)
    GaussianCopulaNVariateProcess(models...)
    GaussianCopulaNVariateProcess(model_1, model_2, ρ₁₂)

Multivariate model joining single name `models` through a Gaussian copula with correlation matrix `ρ`
(the identity when omitted; for two models the scalar correlation `ρ₁₂` can be passed instead).
Each component keeps its own marginal distribution.

# Arguments
- `ρ`: symmetric positive definite correlation matrix, one row per model.
- `models`: the single name models of the components.
"""
struct GaussianCopulaNVariateProcess{num3 <: Number, numtype <: Number} <: NDimensionalMonteCarloProcess{numtype}
    rho::Matrix{num3}
    models::Tuple{Vararg{BaseProcess}}
    function GaussianCopulaNVariateProcess(rho::Matrix{num3}, models::Tuple) where {num3 <: Number}
        check_correlation(rho, models)
        return new{num3, price_type(models..., num3)}(rho, models)
    end
    function GaussianCopulaNVariateProcess(rho::Matrix{num3}, models::BaseProcess...) where {num3 <: Number}
        check_correlation(rho, models)
        return new{num3, price_type(models..., num3)}(rho, models)
    end
    function GaussianCopulaNVariateProcess(models::BaseProcess...)
        len_ = length(models)
        return GaussianCopulaNVariateProcess(Matrix{Float64}(I, len_, len_), models...)
    end
    function GaussianCopulaNVariateProcess(model1::BaseProcess, model2::BaseProcess, rho::num3) where {num3 <: Number}
        corr_matrix_ = [1 rho; rho 1]
        return GaussianCopulaNVariateProcess(corr_matrix_, model1, model2)
    end
end

export GaussianCopulaNVariateProcess;

# The correlation matrix of a copula and its models (#85: one @non_differentiable function, see checks.jl)
function check_correlation(rho, models)
    sz = size(rho)
    check_dimension(sz[1] == sz[2], "The correlation matrix must be square")
    check_dimension(length(models) == sz[1], "One model per row of the correlation matrix")
    check_domain(isposdef(rho), rho, "The correlation matrix must be positive definite")
end
ChainRulesCore.@non_differentiable check_correlation(::Any, ::Any)

function simulate!(S_total, mcProcess::GaussianCopulaNVariateProcess, unds::UnderlyingTuple, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    check_maturity(T)

    ####Simulation
    ## Simulate
    len_ = length(mcProcess.models)
    for i = 1:len_
        simulate!(S_total[i], mcProcess.models[i], unds[i], rfCurve, mcBaseData, T)
    end
    rho = mcProcess.rho
    if (isdiag(rho))
        return
    end
    U_joint = Matrix{eltype(rho)}(undef, len_, Nsim)
    for j = 1:Nstep
        gausscopulagen2!(U_joint, rho, mcBaseData)

        for i = 1:len_
            @views tmp_ = S_total[i][:, j+1]
            # quantile! sorts tmp_ and returns a new vector, so it can be written back into tmp_
            @views tmp_ .= Statistics.quantile!(tmp_, U_joint[i, :])
        end
    end

    ## Conclude
    return
end
