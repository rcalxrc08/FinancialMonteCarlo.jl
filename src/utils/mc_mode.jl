#Parallel Modes
abstract type BaseMode end
struct SerialMode <: BaseMode end
abstract type ParallelMode <: BaseMode end

#Numerical Modes, keyholes for modes different than Monte Carlo
abstract type AbstractMethod end
abstract type AbstractMonteCarloMethod <: AbstractMethod end
struct StandardMC <: AbstractMonteCarloMethod end
struct AntitheticMC <: AbstractMonteCarloMethod end
struct SobolMode <: AbstractMonteCarloMethod end
struct PrescribedMC <: AbstractMonteCarloMethod end