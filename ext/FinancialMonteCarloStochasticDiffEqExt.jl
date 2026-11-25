# Triggered by StochasticDiffEq, which has the SDE solvers since DifferentialEquations 8 is OrdinaryDiffEq only (#64):
# DifferentialEquations 7 loads it, with 8 the user loads it. The initial condition is that of the SDE problem: the
# underlying is not read.
module FinancialMonteCarloStochasticDiffEqExt

using StochasticDiffEq: EnsembleThreads, SOSRI, solve
using FinancialMonteCarlo: AbstractUnderlying, MonteCarloConfiguration, MonteCarloDiffEqModel, ZeroRate, check_maturity
import FinancialMonteCarlo: simulate!

function simulate!(X, mcProcess::MonteCarloDiffEqModel, ::AbstractUnderlying, rfCurve::ZeroRate, mcBaseData::MonteCarloConfiguration, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    Dt = T / Nstep
    check_maturity(T)
    diffeqmodel = mcProcess.model
    sol = solve(diffeqmodel, SOSRI(), EnsembleThreads(), trajectories = Nsim, dt = Dt, adaptive = false)
    X .= [path.u[j] for path in sol.u, j = 1:(Nstep+1)]
    f(x) = mcProcess.final_trasform(x)
    @. X = f(X)
    nothing
end

end
