"""
Struct for VG Process

		vgProcess=VarianceGammaProcess(σ::num1,θ::num2,κ::num3) where {num1,num2,num3<: Number}
	
Where:\n
		σ =	volatility of the process.
		θ = variance of volatility.
		κ =	skewness of volatility.
"""
mutable struct VarianceGammaProcess{num <: Number, num1 <: Number, num2 <: Number, nums0 <: Number, numd <: Number}<:InfiniteActivityProcess
	σ::num
	θ::num1
	κ::num2
	underlying::Underlying{nums0,numd}
	function VarianceGammaProcess(σ::num,θ::num1,κ::num2,underlying::Underlying{nums0,numd}) where {num <: Number, num1 <: Number, num2 <: Number, nums0 <: Number, numd <: Number}
        if σ<=0.0
			error("volatility must be positive");
		elseif κ<=0.0
			error("κappa must be positive");
		elseif 1-σ*σ*κ/2.0-θ*κ<0.0
			error("Parameters with unfeasible values")
		else
            return new{num,num1,num2,nums0,numd}(σ,θ,κ,underlying)
        end
    end
end


export VarianceGammaProcess;

function simulate(mcProcess::VarianceGammaProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	θ1=mcProcess.θ;
	κ1=mcProcess.κ;
	if T<=0.0
		error("Final time must be positive");
	end
	
	dt=T/Nstep;
	#-1/p[3]*log(1+u*u*p[1]*p[1]*p[3]/2.0-1im*p[2]*p[3]*u);
	psi1=-1/κ1*log(1-σ*σ*κ1/2.0-θ1*κ1);
	#1-σ*σ*κ1/2.0-θ1*κ1
	drift=r-d-psi1;
	
	gammaRandomVariable=Gamma(dt/κ1,κ1);
		
	X=simulate(SubordinatedBrownianMotion(σ,drift,gammaRandomVariable,Underlying(0.0)),rfCurve,mcBaseData,T);

	S=S0.*exp.(X);
	
	return S;
	
end
