using .DualNumbers
import Base.Float32
Float32(p::Dual{Float64})=dual(Float32(p.value),Float32(p.epsilon))