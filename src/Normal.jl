# Univariate Normal setup, following Richardson & Green (1997).
# This version uses sampler code that is optimized for univariate normal.
module Normal

module NormalModel # submodule for component family definitions
export Theta, Data, likelihood, log_likelihood, prior_sample, prior_sample!, log_prior,
       construct_hyperparameters, update_parameter!, update_hyperparameters!

include("Random.jl")
using .Random

typealias Theta Array{Float64,1}  # theta = [mu, sigma]
typealias Data Float64

# Normal distribution
const constant = 0.5*log(2.0*pi)
normpdf(x) = exp(-0.5*x*x - constant)
normpdf(x,mu,sigma) = normpdf((x-mu)/sigma)/sigma
log_normpdf(x,mu,sigma) = (r=(x-mu)/sigma;  -0.5*r*r - log(sigma) - constant)
# Gamma distribution
gammapdf(x,a,b) = x^(a-1)*exp(-b*x)*(b^a)/gamma(a)
log_gammapdf(x,a,b) = (a-1)*log(x) - b*x + a*log(b) - lgamma(a)

# Likelihood: Normal(x|mu,sigma^2)
likelihood(x,theta) = normpdf(x,theta[1],theta[2])
likelihood(x,theta,c) = normpdf(x,theta[1,c],theta[2,c])
log_likelihood(x,theta) = log_normpdf(x,theta[1],theta[2])

type Hyperparameters
    m::Float64  # prior mean of mu
    s::Float64  # prior stddev of mu
    a::Float64  # prior shape of sigma
    b::Float64  # prior rate of sigma (starting value)
    g::Float64  # (hyper)prior shape of b
    h::Float64  # (hyper)prior rate of b
end

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

# Prior (Base distribution)
sample_mu(H) = H.s*randn()+H.m
sample_sigma(H) = sqrt(Random.inverse_gamma(H.a,H.b))
prior_sample(H) = [sample_mu(H), sample_sigma(H)]
prior_sample!(theta,n,H) = (for i=1:n; theta[1,i] = sample_mu(H); theta[2,i] = sample_sigma(H); end)
prior_sample!(theta,H) = (theta[1] = sample_mu(H); theta[2] = sample_sigma(H))
log_prior(theta,m,s,a,b) = log_normpdf(theta[1],m,s) + log_gammapdf(1/theta[2]^2,a,b) + log(2/theta[2]^3)
log_prior(theta,H) = log_prior(theta,H.m,H.s,H.a,H.b)

# Update parameter for each component
function update_parameter!(theta_a,theta_b,x,z,c,H,active)
    n = 0
    sum_xc = 0.0
    for i = 1:length(z); if z[i]==c; sum_xc += x[i]; n += 1; end; end
    # update mean
    ell = 1/H.s^2
    lambda = 1/(theta_a[2]^2)
    L = lambda*n + ell
    M = (lambda*sum_xc + ell*H.m) / L
    if active; theta_b[1] = randn()/sqrt(L) + M; end
    
    # update standard deviation
    alpha = H.a + 0.5*n
    r = 0.0
    for i = 1:length(z); if z[i]==c; r += (x[i]-theta_b[1])*(x[i]-theta_b[1]); end; end
    beta = H.b + 0.5*r
    # if active; theta_b[2] = 1/sqrt(rand(Gamma(alpha,1/beta))); end
    if active; theta_b[2] = sqrt(Random.inverse_gamma(alpha,beta)); end

    return log_prior(theta_b,M,1/sqrt(L),alpha,beta)
end

# Update hyperparameters
function update_hyperparameters!(H,theta,ca,t)
    lambda = 1.0./(theta[2,ca[1:t]].^2)
    alpha = H.g + H.a*t
    beta = H.h + sum(lambda)
    H.b = rand(Gamma(alpha, 1/beta))
end

end # module NormalModel
using .NormalModel

# Include generic code
include("generic.jl")

# Include core sampler code
include("coreNormal.jl")

end # module Normal







