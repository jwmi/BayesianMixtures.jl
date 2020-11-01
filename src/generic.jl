
struct Options
    mode::String
    model_type::String
    x::Array{Data,1}
    n_total::Int64
    n_keep::Int64
    n_burn::Int64
    verbose::Bool
    use_hyperprior::Bool
    t_max::Int64 
    # MFM options
    gamma::Float64
    log_pk::String
    # DPM options
    alpha_random::Bool
    alpha::Float64
    # Jain-Neal split-merge options
    use_splitmerge::Bool
    n_split::Int64
    n_merge::Int64
    # RJMCMC options
    k_max::Int64
    # Partition distribution values
    a::Float64
    b::Float64
    log_v::Array{Float64,1}
    # Other
    n::Int64
end

struct Result
    options::Options
    t::Array{Int8,1}
    N::Array{Int16,2}
    z::Array{Int8,2}
    theta::Array{Theta,2}
    keepers::Array{Int64,1}
    elapsed_time::Float64
    time_per_step::Float64
end






