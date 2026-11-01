module FinancialMonteCarloDistributedExt

using Distributed
using FinancialMonteCarlo: AbstractMultiProcess
import FinancialMonteCarlo: default_nbatches, run_on_batches

default_nbatches() = Int64(nworkers())

# The comprehension narrows the Vector{Any} of pmap
run_on_batches(f, configs, ::AbstractMultiProcess) = [result for result in pmap(f, configs)]

end
