"""
    GaussianCopulaNVariateLogProcess(ρ, models...)
    GaussianCopulaNVariateLogProcess(models...)
    GaussianCopulaNVariateLogProcess(model_1, model_2, ρ₁₂)

Like [`GaussianCopulaNVariateProcess`](@ref), but the copula is applied to the log returns
``x ↦ S_0 e^x`` of the components.

# Arguments
- `ρ`: symmetric positive definite correlation matrix, one row per model.
- `models`: the single name models of the components.
"""
struct GaussianCopulaNVariateLogProcess{num3 <: Number, numtype <: Number} <: NDimensionalMonteCarloProcess{numtype}
    models::Tuple{Vararg{BaseProcess}}
    rho::Matrix{num3}
    function GaussianCopulaNVariateLogProcess(rho::Matrix{num3}, models::BaseProcess...) where {num3 <: Number}
        sz = size(rho)
        ChainRulesCore.@ignore_derivatives @assert sz[1] == sz[2]
        ChainRulesCore.@ignore_derivatives @assert length(models) == sz[1]
        ChainRulesCore.@ignore_derivatives @assert isposdef(rho)
        zero_typed = predict_output_type_zero(models...) + zero(num3)
        return new{num3, typeof(zero_typed)}(models, rho)
    end
    function GaussianCopulaNVariateLogProcess(models::BaseProcess...)
        len_ = length(models)
        return GaussianCopulaNVariateLogProcess(Matrix{Float64}(I, len_, len_), models...)
    end
    function GaussianCopulaNVariateLogProcess(model1::BaseProcess, model2::BaseProcess, rho::num3) where {num3 <: Number}
        corr_matrix_ = [1.0 rho; rho 1.0]
        ChainRulesCore.@ignore_derivatives @assert isposdef(corr_matrix_)
        return GaussianCopulaNVariateLogProcess(corr_matrix_, model1, model2)
    end
end

export GaussianCopulaNVariateLogProcess;

function simulate!(S_Total, mcProcess::GaussianCopulaNVariateLogProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    ChainRulesCore.@ignore_derivatives @assert T > 0.0

    ####Simulation
    ## Simulate
    len_ = length(mcProcess.models)

    @inbounds for i = 1:len_
        @views simulate!(S_Total[i], mcProcess.models[i], rfCurve, mcBaseData, T)
    end

    rho = mcProcess.rho
    if (isdiag(rho))
        return
    end
    @inbounds for i = 1:len_
        @. S_Total[i] = log(S_Total[i] / mcProcess.models[i].underlying.S0)
    end
    U_joint = Matrix{eltype(rho)}(undef, len_, Nsim)
    @inbounds for j = 1:Nstep
        gausscopulagen2!(U_joint, rho, mcBaseData)
        for i = 1:len_
            @views tmp_ = S_Total[i][:, j+1]
            sort!(tmp_)
            @views tmp_ .= (mcProcess.models[i].underlying.S0) .* exp.(Statistics.quantile(tmp_, U_joint[i, :]; sorted = true))
        end
    end

    ## Conclude
    return
end
