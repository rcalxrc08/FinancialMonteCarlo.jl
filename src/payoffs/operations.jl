import Base.+;
import Base.*;
import Base.-;

#Strategies Implementation

function +(x::FinancialMonteCarlo.AbstractPayoff,y::FinancialMonteCarlo.AbstractPayoff)
	out=Dict{FinancialMonteCarlo.AbstractPayoff,Number}( x => 1.0 );
	if haskey(out,y)
		out[x]=2.0;
	else
		out[y]=1.0;
	end
	return out;
end

function +(x::Dict{FinancialMonteCarlo.AbstractPayoff,Number},y::FinancialMonteCarlo.AbstractPayoff)
	out=copy(x);
	if haskey(out,y)
		out[y]+=1.0;
	else
		out[y]=1.0;
	end
	return out;
end

function +(x::Dict{FinancialMonteCarlo.AbstractPayoff,Number},y::Dict{FinancialMonteCarlo.AbstractPayoff,Number})
	out=copy(x);
	y_keys=keys(y);
	for y_key in y_keys
		if haskey(out,y_key)
			out[y_key]+=y[y_key];
		else
			out[y_key]=y[y_key];
		end
	end
	return out;
end

function *(x::Dict{FinancialMonteCarlo.AbstractPayoff,Number},y::Number)
	out=copy(x);
	for a in keys(out)
		out[a]*=y;
	end
	return out;
end

function *(x::FinancialMonteCarlo.AbstractPayoff,y::Number)
	return Dict{FinancialMonteCarlo.AbstractPayoff,Number}( x => y );
end
*(y::Number,out::FinancialMonteCarlo.AbstractPayoff)=out*y;
*(y::Number,out::Dict{FinancialMonteCarlo.AbstractPayoff,Number})=out*y;
+(y::FinancialMonteCarlo.AbstractPayoff,out::Dict{FinancialMonteCarlo.AbstractPayoff,Number})= out+y;
-(x::FinancialMonteCarlo.AbstractPayoff)=return -1*x;
-(x::Dict{FinancialMonteCarlo.AbstractPayoff,Number})=return -1*x;
-(y::FinancialMonteCarlo.AbstractPayoff,out::FinancialMonteCarlo.AbstractPayoff)= y+(-1*out);
-(y::FinancialMonteCarlo.AbstractPayoff,out::Dict{FinancialMonteCarlo.AbstractPayoff,Number})= y+(-1*out);
-(y::Dict{FinancialMonteCarlo.AbstractPayoff,Number},out::FinancialMonteCarlo.AbstractPayoff)= y+(-1*out);
-(y::Dict{FinancialMonteCarlo.AbstractPayoff,Number},out::Dict{FinancialMonteCarlo.AbstractPayoff,Number})= y+(-1*out);

#+(x::FinancialMonteCarlo.AbstractPayoff,y::Tuple{FinancialMonteCarlo.AbstractPayoff,FinancialMonteCarlo.AbstractPayoff,typeof(+)})=(x,y,+)

#+(x::Dict{FinancialMonteCarlo.AbstractPayoff,Number},FinancialMonteCarlo.AbstractPayoff)


import Base.Multimedia.display;

function display(p::Dict{FinancialMonteCarlo.AbstractPayoff,Number})
	keys_=keys(p);
	flag=0;
	for key_ in keys_
		val_=p[key_]
		if typeof(val_) <: Real 
			if(flag!=0)
				val_ > 0.0 ? print(+) : print(-);
			end
			if(abs(val_)!=1.0)
				flag!=0 ? print(abs(val_)) : print(val_);
				print(*);
				print(key_);
			else
				print(key_);
			end
		else
			if(flag!=0)
				print(+)
			end
			print("(");
			print(val_);
			print(")");
			print(*);
			print(key_);
		end
		flag+=1;
	end
	println("")
end