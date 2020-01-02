function oper(+,adict)
	@eval import Base.$+
	@eval function $+(r::$adict{num1,num2},d::$adict{num1,num2}) where { num1 <: Number, num2 <: Number }
		T_r=collect(keys_(r));
		T_d=collect(keys_(d));
		if(length(T_d)>length(T_r))
			return $+(d,r)
		end
		d_complete=complete_(T_r,d);
		out=$adict{num1,num2}();
		for t in T_r
			insert!(out,t,  $+(r[t], d_complete[t]))
		end

		return out
	end
	@eval function $+(r::$adict{num1,num2},d::num3) where { num1 <: Number, num2 <: Number , num3 <: Number }
		return $+.(r,d)
	end
end
function oper2(+,adict)
	@eval import Base.$+
	@eval function $+(d::num3,r::$adict{num1,num2}) where { num1 <: Number, num2 <: Number , num3 <: Number }
		return $+(r,d)
	end
end

import Base.-
function -(d::num3,r::Dictionary{num1,num2}) where { num1 <: Number, num2 <: Number , num3 <: Number }
	return +(-1*r,d)
end
function -(d::num3,r::HashDictionary{num1,num2}) where { num1 <: Number, num2 <: Number , num3 <: Number }
	return +(-1*r,d)
end

oper(Symbol(+),FinMCDict)
oper(Symbol(-),FinMCDict)
oper(Symbol(/),FinMCDict)
oper(Symbol(*),FinMCDict)

oper2(Symbol(*),FinMCDict)
oper2(Symbol(+),FinMCDict)


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