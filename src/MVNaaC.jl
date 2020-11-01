# Axis-aligned multivariate normal (i.e., independent entries, i.e., diagonal covariance matrix) using conjugate prior.
module MVNaaC

module MVNaaCmodel # submodule for component family definitions
export Theta, Data, log_marginal, new_theta, Theta_clear!, Theta_adjoin!, Theta_remove!,
       Hyperparameters, construct_hyperparameters, update_hyperparameters!

using Statistics
using SpecialFunctions
lgamma_(x) = logabsgamma(x)[1]

const Data = Array{Float64,1}

mutable struct Theta
    n::Int64                 # number of data points assigned to this cluster
    sum_x::Array{Float64,1}  # sum of the data points x assigned to this cluster
    sum_xx::Array{Float64,1} # sum of x.*x for the data points assigned to this cluster
    Theta(d) = (p=new(); p.n=0; p.sum_x=zeros(d); p.sum_xx=zeros(d); p)
end
new_theta(H) = Theta(H.d)
Theta_clear!(p) = (p.sum_x[:] .= 0.; p.sum_xx[:] .= 0.; p.n = 0)
Theta_adjoin!(p,x) = (for i=1:length(x); p.sum_x[i] += x[i]; p.sum_xx[i] += x[i]*x[i]; end; p.n += 1)
Theta_remove!(p,x) = (for i=1:length(x); p.sum_x[i] -= x[i]; p.sum_xx[i] -= x[i]*x[i]; end; p.n -= 1)

# In each dimension independently,
# X_1,...,X_n ~ Normal(mu,1/lambda) with Normal(mu|m,1/(c*lambda))Gamma(lambda|a,b) prior on mean=mu, precision=lambda.
function log_marginal(p,H)
    n = p.n
    LB = 0.0
    for i=1:H.d; LB += log(H.b + 0.5*p.sum_xx[i] - 0.5*p.sum_x[i]*p.sum_x[i]/n + 0.5*H.c*n*(p.sum_x[i]/n - H.m)^2/(H.c+n)); end
    return H.d*(H.constant - 0.5*n*log(2*pi) - 0.5*log(H.c+n) + H.log_Ga[n]) - (H.a+0.5*n)*LB
end

function log_marginal(x,p,H)
    Theta_adjoin!(p,x)
    result = log_marginal(p,H)
    Theta_remove!(p,x)
    return result
end

mutable struct Hyperparameters
    d::Int64    # dimension
    m::Float64  # prior mean of mu's
    c::Float64  # prior precision multiplier for mu's
    a::Float64  # prior shape of lambda's
    b::Float64  # prior rate of lambda's
    constant::Float64
    log_Ga::Array{Float64,1}
end

function construct_hyperparameters(options)
    x = options.x
    n = length(x)
    d = length(x[1])
    mu = mean(x)
    v = mean([xi.*xi for xi in x]) .- mu.*mu  # sample variance
    @assert(all(abs.(mu) .< 1e-10) && all(abs.(v .- 1.0) .< 1e-10), "Data must be normalized to zero mean, unit variance.")
    m = 0.0
    c = 1.0
    a = 1.0
    b = 1.0
    log_Ga = lgamma_.(a .+ 0.5*(1:n+1))
    constant = 0.5*log(c) + a*log(b) - lgamma_(a)
    return Hyperparameters(d,m,c,a,b,constant,log_Ga)
end

function update_hyperparameters!(H,theta,list,t,x,z)
    #error("update_hyperparameters is not yet implemented.")
end

end # module MVNaaCmodel
using .MVNaaCmodel

# Include generic code
include("generic.jl")

# Include core sampler code
include("coreConjugate.jl")

end # module MVNaaC



