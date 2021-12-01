############### Display Function

import Base.Multimedia.display;
import Base.Multimedia.print;

function print(p::Union{AbstractMonteCarloProcess, AbstractPayoff})
    fldnames = collect(fieldnames(typeof(p)))
    Base.print(typeof(p), "(")
    if (length(fldnames) > 0)
        Base.print(fldnames[1], "=", getfield(p, fldnames[1]))
        popfirst!(fldnames)
        for name in fldnames
            Base.print(",", name, "=", getfield(p, name))
        end
    end
    print(")")
end

function display(p::Union{AbstractMonteCarloProcess, AbstractPayoff})
    fldnames = collect(fieldnames(typeof(p)))
    Base.print(typeof(p), "(")
    if (length(fldnames) > 0)
        Base.print(fldnames[1], "=", getfield(p, fldnames[1]))
        popfirst!(fldnames)
        for name in fldnames
            Base.print(",", name, "=", getfield(p, name))
        end
    end
    println(")")
end

# Agnostic getters and setters, useful for testing and calibration

function get_parameters(model::XProcess) where {XProcess <: AbstractPayoff}
    fields_ = fieldnames(XProcess)
    N1 = length(fields_)
    param_ = []
    for i = 1:N1
        tmp_ = getproperty(model, fields_[i])
        append!(param_, tmp_)
    end
    if (N1 > 0)
        typecheck_ = typeof(sum(param_) + prod(param_))
        param_ = convert(Array{typecheck_}, param_)
    end

    return param_
end

function get_parameters(model::XProcess) where {XProcess <: BaseProcess}
    fields_ = fieldnames(XProcess)
    N1 = length(fields_)
    param_ = []
    for i = 1:N1
        tmp_ = getproperty(model, fields_[i])
        typeof(tmp_) <: AbstractUnderlying ? continue : nothing
        append!(param_, tmp_)
    end
    typecheck_ = typeof(sum(param_) + prod(param_))
    param_ = convert(Array{typecheck_}, param_)

    return param_
end

function get_parameters(model::XProcess) where {XProcess <: VectorialMonteCarloProcess}
    fields_ = fieldnames(XProcess)
    N1 = length(fields_)
    param_ = []
    for i = 1:N1
        tmp_ = getproperty(model, fields_[i])
        typeof(tmp_) <: Tuple ? continue : nothing
        append!(param_, tmp_)
    end
    typecheck_ = typeof(sum(param_) + prod(param_))
    param_ = convert(Array{typecheck_}, param_)

    return param_
end

function set_parameters!(model::XProcess, param::Array{num}) where {XProcess <: BaseProcess, num <: Number}
    fields_ = fieldnames(XProcess)
    N1 = length(fields_) - 1
    @assert N1 == length(param) "Check number of parameters of the model"

    for (symb, nval) in zip(fields_, param)
        if (symb != :underlying)
            setproperty!(model, symb, nval)
        end
    end
end

function set_parameters!(model::LogNormalMixture, param::Array{num}) where {num <: Number}
    N1 = div(length(param) + 1, 2)
    @assert !(N1 * 2 != length(param) + 1 || N1 < 2) "Check number of parameters of the model"
    eta = param[1:N1]
    lam = param[(N1+1):end]
    fields_ = fieldnames(LogNormalMixture)
    setproperty!(model, fields_[1], eta)
    setproperty!(model, fields_[2], lam)
    return
end

function set_parameters!(model::ShiftedLogNormalMixture, param::Array{num}) where {num <: Number}
    N1 = div(length(param), 2)
    @assert !(N1 * 2 != length(param) || N1 < 2) "Check number of parameters of the model"
    eta = param[1:N1]
    lam = param[(N1+1):(2*N1-1)]
    alfa = param[end]
    fields_ = fieldnames(ShiftedLogNormalMixture)
    setproperty!(model, fields_[1], eta)
    setproperty!(model, fields_[2], lam)
    setproperty!(model, fields_[3], alfa)
    return
end
