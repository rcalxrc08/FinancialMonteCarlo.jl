# Type Hierarchy

All models, payoffs and modes provided in this package are organized into a type hierarchy described as follows.

## Models

Remember that each model must have a `simulate` method implemented that specifies the way that the model behaves.
From that abstract class the followings are derived and decomposed as follows:

`AbstractMonteCarloProcess` is used to represent any model that is implemented completely inside this package:
```julia
abstract type AbstractMonteCarloProcess <: AbstractMonteCarloProblem end
```
`ItoProcess` is used to represent any Ito Process:
```julia
abstract type ItoProcess<:AbstractMonteCarloProcess end
```

`LevyProcess` is used to represent any Levy Process:
```julia
abstract type LevyProcess<:AbstractMonteCarloProcess end
```

`FiniteActivityProcess` is used to represent any Finite Activity Levy Process:
```julia
abstract type FiniteActivityProcess<:LevyProcess end
```

`InfiniteActivityProcess` is used to represent any Infinite Activity Levy Process:
```julia
abstract type InfiniteActivityProcess<:LevyProcess end
```

Each concrete model that is implemented inside this package inherits from one of these above models.
Remember that the model type is used to carry the parameters of the model and the dispatch when it's time to simulate, i.e., 
there is just one function `simulate` that is overloaded with respect to the model type.

## Payoffs

The type hierarchy has a root type called `AbstractPayoff`, that type is used whenever it is not necessary to know the model type, like in `pricer`,`var` and `ci`.

Each payoff must have a `payoff` method implemented that specifies the way that the payoff behaves.
From that abstract class the followings are derived and decomposed as follows:

```julia
abstract type AbstractPayoff end
```
```julia
abstract type EuropeanPayoff<:AbstractPayoff end
```

```julia
abstract type AmericanPayoff<:AbstractPayoff end
```

```julia
abstract type BarrierPayoff<:AbstractPayoff end
```

```julia
abstract type AsianPayoff<:AbstractPayoff end
```

Each concrete payoff that is implemented inside this package inherits from one of these above payoffs.
Remember that the payoff type is used to carry the parameters of the model and the dispatch when it's time to price, i.e., 
there is just one function `payoff` that is overloaded with respect to the payoff type.