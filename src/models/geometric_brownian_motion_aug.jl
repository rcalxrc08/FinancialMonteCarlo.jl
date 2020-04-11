"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct GeometricBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number,  num2 <: Number, numtype <: Number} <: ItoProcess{numtype}
	σ::num
	μ::Curve{num1,num4}
	x0::num2
	function GeometricBrownianMotionVec(σ::num,μ::Curve{num1,num4},x0::num2) where {num <: Number , num1 <: Number, num4 <: Number, num2 <: Number}
        if σ <= 0
            error("Volatility must be positive")
        else
			zero_typed=zero(num)+zero(num1)+zero(num4)+zero(num2);
            return new{num,num1,num4,num2,typeof(zero_typed)}(σ,μ,x0)
        end
    end
end

function GeometricBrownianMotion(σ::num,μ::Curve{num1,num4},x0::num2) where {num <: Number, num1 <: Number, num4 <: Number, num2 <: Number}
	return GeometricBrownianMotionVec(σ,μ,x0)
end

export GeometricBrownianMotionVec;

function simulate!(X,mcProcess::GeometricBrownianMotionVec,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	@assert T>0
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-(σ_gbm^2/2);
	simulate!(X,BrownianMotion(σ_gbm,μ_bm),mcBaseData,T)
	S0=mcProcess.x0;
	f(x)=S0*exp(x);
	broadcast!(f,X,X)
	nothing;
end