using IterTools
#My implementation of lexicographic less function, not sure if it is the most performant but it works
function lex_less(x, y)
    if (length(x) != length(y))
        return length(y) > length(x)
    else
        x_s = sort(x)
        y_s = sort(y)
        for i = 1:length(x)
            if (x_s[i] != y_s[i])
                return y_s[i] > x_s[i]
            end
        end
    end
end

#This function solves the vectorial problem:
#Multi dimensional models have no clue of the name of the underlying process and no clue of the name of themselves,
#A pricer implemented for a multidimensional model with a single name option makes little sense since there his no
#way to map the option to the right underlying.
#In order to solve this problem the pricer is implemented in a strange way:
#A vector of integer is computed in the following way:
#Assuming 3 underlyings
#1st element ---> option referred to model[1]
#2nd element ---> option referred to model[2]
#3rd element ---> option referred to model[3]
#4th element ---> option referred to (model[1],model[2])
#5th element ---> option referred to (model[1],model[3])
#6th element ---> option referred to (model[2],model[3])
#7th element ---> option referred to (model[1],model[2],model[3])
function compute_indices(x::Integer)
    X = collect(1:x)
    #Subsets will return the empty set as well
    f(y) = filter(x -> !(length(x) == 0), y)
    k = subsets(X) |> collect |> f
    out = sort(k, lt = (x, y) -> lex_less(x, y))

    return out
end

#Following function returns the names (duplicated as well) of the underlying names that must simulate.
#An underlying to be simulated must have an option in the portfolio that is referred to it (directly or in a basket).
function complete_2(payoffs_underlyings, models_underlyings)
    out = String[]
    for x_ in payoffs_underlyings
        for y_ in models_underlyings
            if (any(z -> z == x_, split(y_, "_")) || y_ == x_)
                push!(out, y_)
            end
        end
    end
    return out
end
