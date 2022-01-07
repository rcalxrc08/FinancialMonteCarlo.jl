# # Getting Started
#
# ## Installation
# The FinancialMonteCarlo package is available through the Julia package system 
# by running `]add FinancialMonteCarlo`.
# Throughout, we assume that you have installed the package.

# ## Basic syntax
# The basic syntax for FinancialMonteCarlo is simple. Here is an example of pricing a European Option with a Black Scholes model:

using FinancialMonteCarlo
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;
# Define FinancialMonteCarlo Parameters:
Nsim = 10000;
Nstep = 30;
# Define Model Parameters.
σ = 0.2;
# Build the MonteCarlo configuration and the zero rate.
mcConfig = MonteCarloConfiguration(Nsim, Nstep);
rfCurve = ZeroRate(r);
# Define The Option
EuOption = EuropeanOption(T, K)
# Define the Model of the Underlying
Model = BlackScholesProcess(σ, Underlying(S0, d));
# Call the pricer metric
@show EuPrice = pricer(Model, rfCurve, mcConfig, EuOption);
# ## Using Other Models and Options
# The package contains a large number of models of three main types:
# * `Ito Process`
# * `Jump-Diffusion Process`
# * `Infinity Activity Process`

# And the following types of options:

# * `European Options`
# * `Asian Options`
# * `Barrier Options`
# * `American Options`

# The usage is analogous to the above example.