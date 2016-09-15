# Univariate Normal setup, following Richardson & Green (1997).
# This version uses the general version of the sampler.
module NormalGeneral

module NormalModel
export Theta, Data, log_likelihood, log_prior, prior_sample!, new_theta, Theta_clear!, Theta_adjoin!, Theta_remove!,
       Hyperparameters, construct_hyperparameters, update_hyperparameters!, update_parameter!

include("Gamma.jl")
using .Gamma

typealias Data Float64

type Normal_params
    mu::Float64
    sigma::Float64
    n::Int64        # number of data points assigned to this cluster
    sum_x::Float64  # sum of the data points x assigned to this cluster
    sum_xx::Float64 # sum of x^2 for the data points assigned to this cluster
    Normal_params(mu,sigma) = (p=new(); p.mu=mu; p.sigma=sigma; p.n=0; p.sum_x=0.0; p.sum_xx=0.0; p)
end

Normal_logpdf(x,mu,sigma) = (r=(x-mu)/sigma;  -0.5*r*r - 0.5*log(2.0*pi) - log(sigma))
Normal_clear!(p) = (p.sum_x = 0.; p.sum_xx = 0.; p.n = 0)
Normal_adjoin!(p,x) = (p.sum_x += x; p.sum_xx += x*x; p.n += 1)
Normal_remove!(p,x) = (p.sum_x -= x; p.sum_xx -= x*x; p.n -= 1)

type Hyperparameters
    m::Float64  # prior mean of mu
    s::Float64  # prior stddev of mu
    a::Float64  # prior shape of sigma
    b::Float64  # prior rate of sigma (starting value)
    g::Float64  # (hyper)prior shape of b
    h::Float64  # (hyper)prior rate of b
end

Gamma_logpdf(x,a,b) = (a-1)*log(x) - b*x + a*log(b) - lgamma(a)
# log density of SqrtInvGamma (distn of sigma when 1/sigma^2 ~ Gamma(a,b))
SqrtInvGamma_logpdf(x,a,b) = Gamma_logpdf(1/x^2,a,b) + log(2/x^3)

typealias Theta Normal_params
log_likelihood(x,p) = Normal_logpdf(x,p.mu,p.sigma)
# prior: Normal(mu|H.m,H.s) * SqrtInvGamma(sigma|H.a,H.b)
log_prior(p,m,s,a,b) = Normal_logpdf(p.mu,m,s) + SqrtInvGamma_logpdf(p.sigma,a,b)
log_prior(p,H) = log_prior(p,H.m,H.s,H.a,H.b)
prior_sample!(p,H) = (p.mu = randn()*H.s + H.m; p.sigma = 1/sqrt(Gamma.gamma(H.a,H.b)))
new_theta(H) = Normal_params(0.0,1.0)
Theta_clear!,Theta_adjoin!,Theta_remove! = Normal_clear!,Normal_adjoin!,Normal_remove!


function construct_hyperparameters(options)
    x = options.x
    # Use values in Green & Richardson (1997) to enable comparison.
    m = (minimum(x) + maximum(x))/2
    R = maximum(x) - minimum(x)
    s = R
    a = 2.0
    b = 1.0  # b will be updated in Gibbs sampling
    g = 0.2
    h = 10/R^2
    return Hyperparameters(m,s,a,b,g,h)
end


function update_parameter!(theta_a,theta_b,H,active,density)
    n,sum_x,sum_xx = theta_b.n,theta_b.sum_x,theta_b.sum_xx
    # update mean
    ell = 1/H.s^2
    lambda = 1/(theta_a.sigma^2)
    L = lambda*n + ell
    M = (lambda*sum_x + ell*H.m) / L
    if active; theta_b.mu = randn()/sqrt(L) + M; end
    
    # update standard deviation
    alpha = H.a + 0.5*n
    beta = H.b + 0.5*(sum_xx - 2*theta_b.mu*sum_x + n*theta_b.mu^2)
    if active; theta_b.sigma = 1/sqrt(Gamma.gamma(alpha,beta)); end

    return (density? log_prior(theta_b,M,1/sqrt(L),alpha,beta) : NaN)
end

function update_hyperparameters!(H,theta,list,t,x,z)
    sum_lambdas = 0.0
    for k = 1:t; sum_lambdas += 1/theta[list[k]].sigma^2; end
    alpha = H.g + H.a*t
    beta = H.h + sum_lambdas
    H.b = Gamma.gamma(alpha,beta)
end

end # module NormalModel
using .NormalModel

# Include generic code
include("generic.jl")

# Include core sampler code
include("coreNonconjugate.jl")

end # module NormalGeneral






