function oper(+)
	@eval import Base.$+
	@eval function $+(r::FinMCDict{num1,num2},d::FinMCDict{num1,num3}) where { num1 <: Number, num2 <: Number, num3 <: Number }
		T_r=collect(keys_(r));
		T_d=collect(keys_(d));
		if(length(T_d)>length(T_r))
			return $+(d,r)
		end
		d_complete=complete_(T_r,d);
		#key_type=typeof(zero(num1));
		val_type=typeof(zero(num2)+zero(num3));
		out=FinMCDict{num1,val_type}();
		for t in T_r
			insert!(out,t,  $+(r[t], d_complete[t]))
		end

		return out
	end
	@eval function $+(r::FinMCDict{num1,num2},d::num3) where { num1 <: Number, num2 <: Number , num3 <: Number }
		return $+.(r,d)
	end
end
function oper2(+)
	@eval import Base.$+
	@eval function $+(d::num3,r::FinMCDict{num1,num2}) where { num1 <: Number, num2 <: Number , num3 <: Number }
		return $+(r,d)
	end
end

import Base.-
function -(d::num3,r::FinMCDict{num1,num2}) where { num1 <: Number, num2 <: Number , num3 <: Number }
	return +(-1*r,d)
end

oper(Symbol(+))
oper(Symbol(-))
oper(Symbol(/))
oper(Symbol(*))

oper2(Symbol(*))
oper2(Symbol(+))


keys_(x::Dictionary)=x.indices		  
keys_(x)=keys(x)

values_(x::Dictionary)=x.values		  
values_(x)=values(x)

function complete_(T,d::FinMCDict{num1,num2}) where {num1,num2}
	T_d=collect(keys_(d));
	idx_=sortperm(T_d);
	T_d=T_d[idx_]
	d_val=collect(values_(d));
	d_val=d_val[idx_]
	out=FinMCDict{num1,num2}();
	itp=LinearInterpolation(T_d,d_val);
	for t in T
		insert!(out,t,itp(t))
	end
	return out;
end