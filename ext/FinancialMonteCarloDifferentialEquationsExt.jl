module FinancialMonteCarloDifferentialEquationsExt

using ChainRulesCore
using DifferentialEquations: EnsembleThreads, SOSRI, solve
using FinancialMonteCarlo: MonteCarloConfiguration, MonteCarloDiffEqModel, ZeroRate, check_domain
import FinancialMonteCarlo: simulate!

function simulate!(X, mcProcess::MonteCarloDiffEqModel, rfCurve::ZeroRate, mcBaseData::MonteCarloConfiguration, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    Dt = T / Nstep
    ChainRulesCore.@ignore_derivatives check_domain(T > 0, T, "Time to Maturity must be positive")
    diffeqmodel = mcProcess.model
    sol = solve(diffeqmodel, SOSRI(), EnsembleThreads(), trajectories = Nsim, dt = Dt, adaptive = false)
    X .= [path.u[j] for path in sol.u, j = 1:(Nstep+1)]
    f(x) = mcProcess.final_trasform(x)
    @. X = f(X)
    nothing
end

end
