import Base.+;

#Strategies Implementation

const Portfolio = Dict{String, Position};

#Link an underlying name to a Position
function →(x::String, y::Position)
    out = Portfolio(x => y)

    return out
end

#Link an underlying name to a SingleNamePayoff
function →(x::String, y::SingleNamePayoff)
    sep = findfirst("_", x)
    @assert isnothing(sep) "NO UNDERSCORE ALLOWED IN SINGLE NAME OPTIONS"
    return x → (1.0 * y)
end

#Link an underlying name to a BasketPayoff
function →(x::String, y::BasketPayoff)
    sep = findfirst("_", x)[1]
    @assert !isnothing(sep) "Bivariate payoff underlying must follow the format: INDEX1_INDEX2"
    return x → (1.0 * y)
end

#Join Portfolio s
function +(x::Portfolio, y::Portfolio)
    out = copy(x)
    y_keys = keys(y)
    for y_key in y_keys
        if haskey(out, y_key)
            out[y_key] = out[y_key] + y[y_key]
        else
            out[y_key] = y[y_key]
        end
    end
    return out
end

#Overload for display
import Base.Multimedia.display;
function display(p::Portfolio)
    keys_ = keys(p)
    for key_ in keys_
        print(key_)
        print(" => ")
        print(p[key_])
    end
    println("")
end

#Overload for extraction, output needs particular shape.
function extract_(x::String, dict_::Portfolio)
    keys_ = keys(dict_)
    out = Position()
    out_vec = Position[]
    ap_keys = String[]
    for key_ in keys_
        if ((x == key_) || (any(z -> z == key_, split(x, "_"))))
            out = dict_[key_]
            push!(out_vec, out)
            push!(ap_keys, key_)
        end
    end

    if (length(out_vec) == 1) && (x == ap_keys[1])
        return out
    else
        str_v = split(x, "_")
        X = compute_indices(length(str_v))

        #tmp_map=Dict( str_v .=> collect(1:length(str_v)))
        tmp_map_rev = Dict(collect(1:length(str_v)) .=> str_v)

        out__ = Array{Position}(undef, length(X))

        for i_ = 1:length(X)
            idx_ = X[i_]
            str_underlv = [tmp_map_rev[idx_2] for idx_2 in idx_]
            str_underl = join(str_underlv, "_")

            if (haskey(dict_, str_underl))
                out__[i_] = dict_[str_underl]
            end
        end

        return out__
    end
end
