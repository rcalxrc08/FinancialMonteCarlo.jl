# simulate! is in the DifferentialEquations extension.
struct MonteCarloDiffEqModel{absdiffeqmodel} <: ScalarMonteCarloProcess{Float64}
    model::absdiffeqmodel
    final_trasform::Function
    underlying::AbstractUnderlying
    function MonteCarloDiffEqModel(model::absdiffeqmodel, final_trasform::Function, underlying::AbstractUnderlying) where {absdiffeqmodel}
        return new{absdiffeqmodel}(model, final_trasform, underlying)
    end
    function MonteCarloDiffEqModel(model::absdiffeqmodel, underlying::AbstractUnderlying) where {absdiffeqmodel}
        func(x) = identity(x)
        return new{absdiffeqmodel}(model, func, underlying)
    end
end

export MonteCarloDiffEqModel;
