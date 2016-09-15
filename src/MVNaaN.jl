# Axis-aligned multivariate normal (i.e., independent entries, i.e., diagonal covariance matrix)
# using conjugate prior, but not integrating it out. This is just to enable comparison with MVNaaC.
module MVNaaN

module MVNaaNmodel # submodule for component family definitions
export Theta, Data, log_likelihood, log_prior, prior_sample!, new_theta, Theta_clear!, Theta_adjoin!, Theta_remove!,
       Hyperparameters, construct_hyperparameters, update_hyperparameters!, update_parameter!

include("Gamma.jl")
using .Gamma

typealias Data Array{Float64,1}

type Theta
    mu::Array{Float64,1}     # means
    lambda::Array{Float64,1} # precisions
    d::Int64                 # dimension
    n::Int64                 # number of data points assigned to this cluster
    sum_x::Array{Float64,1}  # sum of the data points x assigned to this cluster
    sum_xx::Array{Float64,1} # sum of x.*x for the data points assigned to this cluster
    Theta(mu,lambda) = (p=new(); p.mu=mu; p.lambda=lambda; p.d=length(mu); p.n=0; p.sum_x=zeros(p.d); p.sum_xx=zeros(p.d); p)
end
new_theta(H) = Theta(zeros(H.d),ones(H.d))

const constant = -0.5*log(2*pi)
Normal_logpdf(x,mu,lambda) = -0.5*lambda*(x-mu)*(x-mu) + 0.5*log(lambda) + constant
MVNaaN_logpdf(x,mu,lambda) = (r=0.0; for i=1:length(x); r+=Normal_logpdf(x[i],mu[i],lambda[i]); end; r)
Theta_clear!(p) = (for i=1:p.d; p.sum_x[i] = 0.0; p.sum_xx[i] = 0.0; end; p.n = 0)
Theta_adjoin!(p,x) = (for i=1:p.d; p.sum_x[i] += x[i]; p.sum_xx[i] += x[i]*x[i]; end; p.n += 1)
Theta_remove!(p,x) = (for i=1:p.d; p.sum_x[i] -= x[i]; p.sum_xx[i] -= x[i]*x[i]; end; p.n -= 1)

type Hyperparameters
    d::Int64    # dimension
    m::Float64  # prior mean of mu's
    c::Float64  # prior precision multiplier for mu's
    a::Float64  # prior shape of lambda's
    b::Float64  # prior rate of lambda's
end

Gamma_logpdf(x,a,b) = (a-1)*log(x) - b*x + a*log(b) - lgamma(a)

log_likelihood(x,p) = MVNaaN_logpdf(x,p.mu,p.lambda)
# prior: Normal(mu[i]|H.m,1/(H.c*lambda[i])) * Gamma(lambda[i]|H.a,H.b)
log_prior(p,m,c,a,b) = (r=0.0; for i=1:p.d; r += Normal_logpdf(p.mu[i],m,c*p.lambda[i]) + Gamma_logpdf(p.lambda[i],a,b); end; r)
log_prior(p,H) = log_prior(p,H.m,H.c,H.a,H.b)
prior_sample!(p,H) = (for i=1:H.d; p.lambda[i] = Gamma.gamma(H.a,H.b); p.mu[i] = randn()/sqrt(H.c*p.lambda[i]) + H.m; end)


function construct_hyperparameters(options)
    x = options.x
    d = length(x[1])
    mu = mean(x)
    v = mean([xi.*xi for xi in x]) - mu.*mu  # sample variance
    @assert(all(abs(mu) .< 1e-10) && all(abs(v - 1.0) .< 1e-10), "Data must be normalized to zero mean, unit variance.")
    m = 0.0
    c = 1.0
    a = 1.0
    b = 1.0
    return Hyperparameters(d,m,c,a,b)
end

function update_parameter!(theta_a,theta_b,H,active,density)
    n,sum_x,sum_xx = theta_b.n,theta_b.sum_x,theta_b.sum_xx
    logp = 0.0
    A = H.a + 0.5*(n+1)
    for i = 1:H.d
        # update mean
        L = theta_a.lambda[i]*(n + H.c)
        M = (sum_x[i] + H.c*H.m) / (n + H.c)
        if active; theta_b.mu[i] = randn()/sqrt(L) + M; end
        if density; logp += Normal_logpdf(theta_b.mu[i],M,L); end
        
        # update precision
        B = H.b + 0.5*(sum_xx[i] - 2*theta_b.mu[i]*sum_x[i] + n*theta_b.mu[i]^2) + 0.5*H.c*(theta_b.mu[i] - H.m)^2
        if active; theta_b.lambda[i] = Gamma.gamma(A,B); end
        if density; logp += Gamma_logpdf(theta_b.lambda[i],A,B); end
    end
    return (density? logp : NaN)
end

function update_hyperparameters!(H,theta,list,t,x,z)
    #error("update_hyperparameters is not yet implemented.")
end

end # module MVNaaNmodel
using .MVNaaNmodel

# Include generic code
include("generic.jl")

# Include core sampler code
include("coreNonconjugate.jl")

end # module MVNaaN



