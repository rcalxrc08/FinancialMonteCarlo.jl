# simulate! is in the StochasticDiffEq extension.
struct MonteCarloDiffEqModel{absdiffeqmodel} <: ScalarMonteCarloProcess{Float64}
    model::absdiffeqmodel
    final_trasform::Function
    function MonteCarloDiffEqModel(model::absdiffeqmodel, final_trasform::Function) where {absdiffeqmodel}
        return new{absdiffeqmodel}(model, final_trasform)
    end
    function MonteCarloDiffEqModel(model::absdiffeqmodel) where {absdiffeqmodel}
        func(x) = identity(x)
        return new{absdiffeqmodel}(model, func)
    end
end

export MonteCarloDiffEqModel;
