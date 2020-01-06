using IterTools
function lex_less(x,y)
	if(length(x)!=length(y))
		return length(y) > length(x);
	else
		x_s=sort(x)
		y_s=sort(y)
		for i in 1:length(x)
			if(x_s[i] != y_s[i])
				return y_s[i] > x_s[i];
			end
		end
	end
end

function compute_indices(x::Integer)

	X=collect(1:x);
	f(y)=filter(x->!(length(x)==0),y)
	k=subsets(X)|>collect|>f
	out=sort(k,lt=(x,y)->lex_less(x,y))
	
	return out;
	
end

function complete_2(x,y)
	out=String[];
	for x_ in x
		for y_ in y
			if(any(z->z==x_,split(y_,"_"))||y_==x_)
				push!(out,y_)
			end
		end
	end
	return out;
end