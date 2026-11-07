module FinancialMonteCarloDifferentialEquationsExt

using ChainRulesCore
using DifferentialEquations: EnsembleThreads, SOSRI, solve
using FinancialMonteCarlo: MonteCarloConfiguration, MonteCarloDiffEqModel, ZeroRate
import FinancialMonteCarlo: simulate!

function simulate!(X, mcProcess::MonteCarloDiffEqModel, rfCurve::ZeroRate, mcBaseData::MonteCarloConfiguration, T::numb) where {numb <: Number}
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

end
