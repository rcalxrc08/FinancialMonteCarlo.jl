abstract type AbstractSimResult end

struct SimResult{num_ <: Number, num_2 <: Number, time_array <: AbstractArray{num_}, sim_matrix <: AbstractMatrix{num_2}} <: AbstractSimResult
    times::time_array
    sim_result::sim_matrix
    initial_value::num_2
end
