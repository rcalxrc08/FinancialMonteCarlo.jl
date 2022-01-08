function oper(+)
    @eval import Base.$+
    @eval function $+(r::DictTypeInternal{num1, num2}, d::DictTypeInternal{num1, num3}) where {num1 <: Number, num2 <: Number, num3 <: Number}
        T_r = extract_keys(r)
        T_d = extract_keys(d)
        if (length(T_d) > length(T_r))
            return $+(d, r)
        end
        d_complete = complete_curve(T_r, d)
        val_type = typeof(zero(num2) + zero(num3))
        out = DictTypeInternal{num1, val_type}()
        for t in T_r
            insert!(out, t, $+(r[t], d_complete[t]))
        end

        return out
    end
    @eval function $+(r::DictTypeInternal{num1, num2}, d::num3) where {num1 <: Number, num2 <: Number, num3 <: Number}
        return $+.(r, d)
    end
end
function oper2(+)
    @eval import Base.$+
    @eval function $+(d::num3, r::DictTypeInternal{num1, num2}) where {num1 <: Number, num2 <: Number, num3 <: Number}
        return $+(r, d)
    end
end

import Base.-
function -(d::num3, r::DictTypeInternal{num1, num2}) where {num1 <: Number, num2 <: Number, num3 <: Number}
    return +(-1 * r, d)
end

oper(Symbol(+))
oper(Symbol(-))
oper(Symbol(/))
oper(Symbol(*))

oper2(Symbol(*))
oper2(Symbol(+))

function complete_curve(T::Array{num}, d::DictTypeInternal{num1, num2}) where {num <: Number, num1 <: Number, num2 <: Number}
    T_d = extract_keys(d)
    idx_ = sortperm(T_d)
    T_d = T_d[idx_]
    d_val = extract_values(d)
    d_val = d_val[idx_]
    out = DictTypeInternal{num1, num2}()
    itp = LinearInterpolation(T_d, d_val)
    for t in T
        insert!(out, t, itp(t))
    end
    return out
end
