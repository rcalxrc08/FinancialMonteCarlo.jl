include("../testutils.jl")
using ChainRulesTestUtils

# The reverse rules of the package (#80), against finite differences: integral of a curve, set_parameters!
@testset "ChainRules" begin
    @testset "integral of a curve, x = $x" for x in (0.1, 0.6, 0.7, 1.5, 3.0) # not on a node, where the integral has a kink
        curve = FinancialMonteCarlo.Curve([0.01, 0.03, 0.02, 0.05], [0.2, 0.5, 1.0, 2.0])
        test_rrule(FinancialMonteCarlo.integral, curve, x)
        test_rrule(FinancialMonteCarlo.integral, FinancialMonteCarlo.Curve([0.02], [1.0]), x)
    end
    underlying = Underlying(S0, FinancialMonteCarlo.ImpliedCurve([0.01, 0.012, 0.015], [0.5, 1.0, 2.0]))
    models = (BlackScholesProcess(sigma, underlying), HestonProcess(0.3, 0.2, 0.0, 1.5, -0.5, 0.04, underlying), LogNormalMixture((0.2, 0.3), (0.6,), underlying), ShiftedProcess(LogNormalMixture((0.2, 0.3, 0.25), (0.4, 0.3), underlying), 0.1), ShiftedProcess(BlackScholesProcess(sigma, Underlying(S0, d)), 0.2))
    @testset "set_parameters! $(nameof(typeof(model)))" for model in models
        test_rrule(FinancialMonteCarlo.set_parameters!, model, FinancialMonteCarlo.get_parameters(model) .* 1.01)
    end
end
