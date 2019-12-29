"""
Struct for Kou Process

		kouProcess=KouProcess(σ::num1,λ::num2,p::num3,λ₊::num4,λ₋::num5) where {num1,num2,num3,num4,num5 <: Number}
	
Where:\n
		σ  =	volatility of the process.
		λ  = 	jumps intensity.
		p  =	prob. of having positive jump.
		λ₊ =	positive jump size.
		λ₋ =	negative jump size.
"""
mutable struct KouProcess{num <: Number, num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number, nums0 <: Number, numd <: Number}<:FiniteActivityProcess
	σ::num
	λ::num1
	p::num2
	λ₊::num3
	λ₋::num4
	underlying::Underlying{nums0,numd}
	function KouProcess(σ::num,λ::num1,p::num2,λ₊::num3,λ₋::num4,underlying::Underlying{nums0,numd}) where {num <: Number, num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number, nums0 <: Number, numd <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif λ<=0.0
			error("jump intensity must be positive");
		elseif λ₊<=0.0
			error("positive λ must be positive");
		elseif λ₋<=0.0
			error("negative λ must be positive");
		elseif !(0<=p<=1)
			error("p must be a probability")
        else
            return new{num,num1,num2,num3,num4,nums0,numd}(σ,λ,p,λ₊,λ₋,underlying)
        end
    end
end

export KouProcess;

compute_jump_size(mcProcess::KouProcess,mcBaseData::MonteCarloConfiguration)=rand(mcBaseData.rng)<mcProcess.p ? quantile_exp(mcProcess.λ₊,rand(mcBaseData.rng)) : -quantile_exp(mcProcess.λ₋,rand(mcBaseData.rng))
compute_drift(mcProcess::KouProcess)=-(-mcProcess.σ^2/2-mcProcess.λ*(mcProcess.p/(mcProcess.λ₊-1)-(1-mcProcess.p)/(mcProcess.λ₋+1)))