module FinancialMonteCarloDifferentialEquationsExt

using ChainRulesCore
using DifferentialEquations: EnsembleThreads, SOSRI, solve
using FinancialMonteCarlo: AbstractMonteCarloMethod, BaseMode, MonteCarloConfiguration, MonteCarloDiffEqModel, ZeroRate
import FinancialMonteCarlo: predict_output_type_zero, simulate!

function simulate!(X, mcProcess::MonteCarloDiffEqModel, rfCurve::ZeroRate, mcBaseData::MonteCarloConfiguration{type1, type2, type3, type4}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    Dt = T / Nstep
    ChainRulesCore.@ignore_derivatives @assert T > 0.0
    diffeqmodel = mcProcess.model
    sol = solve(diffeqmodel, SOSRI(), EnsembleThreads(), trajectories = Nsim, dt = Dt, adaptive = false)
    X .= [path.u[j] for path in sol.u, j = 1:(Nstep+1)]
    f(x) = mcProcess.final_trasform(x)
    @. X = f(X)
    nothing
end

function predict_output_type_zero(::absdiffeqmodel) where {absdiffeqmodel <: MonteCarloDiffEqModel}
    return zero(Float64)
end

end
