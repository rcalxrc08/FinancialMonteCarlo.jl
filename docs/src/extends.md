# Create New Models and Payoffs

Whereas this package already provides a large collection of common models/payoffs out of box, there are still occasions where you want to create new models/payoffs (*e.g* your application requires a special kind of models/payoffs, or you want to contribute to this package).

Generally, you don't have to implement every API method listed in the documentation. This package provides a series of generic functions that turn a small number of internal methods into user-end API methods. What you need to do is to implement this small set of internal methods for your models/payoffs.

## Create a Model
To use your model inside the package, you must implement a struct, that is inheriting from one of the abstract class define in the Type section.
Such struct must contain the parameters of your model, then you have to implement a function called `simulate` that is representing how your model simulates.
You can easily copy the signature of the function from one of the implementation in the models folder (beside subordinated brownian motion).

## Create a Payoff
To use your payoff inside the package, you must implement a struct, that is inheriting from one of the abstract class define in the Type section.
Such struct must contain the parameters (strike, time to maturity, barriers, etc) of your payoff, then you have to implement a function called `simulate` that is representing how your payoff simulates.
You can easily copy the signature of the function from one of the implementation in the payoffs folder (beside general american option).

## Finalize the interaction
Once you have implemented your own payoff/model, you can directly use the functionalities provided already by the package, e.g. `pricer`,`var`,`ci` and all of the models/payoffs that have been already implemented.