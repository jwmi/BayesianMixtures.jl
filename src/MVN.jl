# Multivariate normal setup using non-conjugate prior.
module MVN

module MVNmodel
export Theta, Data, log_likelihood, log_prior, prior_sample!, new_theta, Theta_clear!, Theta_adjoin!, Theta_remove!,
       Hyperparameters, construct_hyperparameters, update_hyperparameters!, update_parameter!, mixrnd, mixture_density
       
include("Lower.jl")
using .Lower

const NU_SIGMA_PROP = 0.1  # Metropolis proposal std dev for H.R.nu updates
    
typealias Data Array{Float64,1}

# Multivariate normal distribution
type MVN_params
    # read/write
    m::Array{Float64,1} # mean
    # read/restricted-write (must notify after writing)
    L::Array{Float64,2} # lower triangular matrix such that L*L' = R
    # no access (use getter to read when necessary)
    _R::Array{Float64,2} # precision matrix (inverse covariance)
    _R_valid::Bool # flag to indicate whether R has been computed or not
    # read-only
    logdetR::Float64 # log of the determinant of the precision matrix
    d::Int64 # dimension
    n::Int64 # number of data points assigned to this cluster
    sum_x::Array{Float64,1} # sum of the data points x assigned to this cluster
    sum_xx::Array{Float64,2} # sum of x*x' for the data points assigned to this cluster
    function MVN_params(m,R)
        p = new(); p.m = copy(m); p._R = copy(R); p._R_valid = true; p.d = d = length(m)
        p.L = zeros(d,d); Lower.Cholesky!(p.L,p._R,d)
        p.logdetR = Lower.logdetsq(p.L,d)
        p.n = 0
        p.sum_x = zeros(d)
        p.sum_xx = zeros(d,d)
        return p
    end
end
# Interface functions for MVN_params
MVN_logpdf(x,p) = -0.5*Lower.quadratic(x, p.m, p.L, p.d) - 0.5*p.d*log(2*pi) + 0.5*p.logdetR
MVN_sample!(x,p) = Lower.sample_Normal!(x, p.m, p.L, p.d)
MVN_get_R!(p) = (if !p._R_valid; Lower.multiplyMNt!(p._R,p.L,p.L,p.d); p._R_valid=true; end; p._R)
MVN_notify_L!(p) = (p.logdetR = Lower.logdetsq(p.L,p.d); p._R_valid = false) # call this after modifying L
MVN_clear!(p) = (fill!(p.sum_x,0.); fill!(p.sum_xx,0.); p.n=0)
MVN_adjoin!(p,x) = (for i=1:p.d; p.sum_x[i]+=x[i]; for j=1:p.d; p.sum_xx[i,j]+=x[i]*x[j]; end; end; p.n+=1)
MVN_remove!(p,x) = (for i=1:p.d; p.sum_x[i]-=x[i]; for j=1:p.d; p.sum_xx[i,j]-=x[i]*x[j]; end; end; p.n-=1)



# Wishart distribution
type Wishart_params
    # read/restricted-write (must notify after writing)
    W::Array{Float64,2} # inverse of the mean matrix, i.e., W = inv(nu*V)
    # read-only
    nu::Float64 # degrees of freedom (can modify using set_nu)
    M::Array{Float64,2} # lower triangular matrix such that M*M' = inv(nu*W) = V = the scale matrix
    c::Float64 # log normalization constant
    d::Int64 # dimension
    n::Int64 # number of observations
    sum_R::Array{Float64,2} # sum of the R's
    sum_logdetR::Float64 # sum of the log(det(R))'s
    A::Array{Float64,2} # temporary variable
    function Wishart_params(nu,W) # W is the inverse of the mean matrix (d-by-d), nu > d-1 is the degrees of freedom
        p = new(); p.nu = nu; p.W = copy(W); p.d = d = size(p.W,1)
        p.M = zeros(d,d); p.A = zeros(d,d)
        Lower.Cholesky_inverse!(p.M,p.A,p.W,d); Lower.scale!(p.M,sqrt(1/p.nu),d)
        p.c = Wishart_constant(p)
        p.n = 0
        p.sum_R = zeros(d,d)
        p.sum_logdetR = 0.
        return p
    end
end
# Helper functions
matrix_dot(A,B,d) = (r=0.; for i=1:d*d; r += A[i]*B[i]; end; r)
Wishart_constant(p) = (d=p.d; nu=p.nu; r=0.; for j=1:d; r += lgamma(0.5*(nu+1-j)); end;
    0.5*nu*Lower.logdetsq(p.M,d) + 0.5*nu*d*log(2) + 0.25*d*(d-1)*log(pi) + r)
# Interface functions
# Sample L (lower triangular) such that R = L*L' ~ Wishart(Scale=M*M', DOF=nu)
Wishart_sample!(L,p) = Lower.sample!(L, p.M, p.nu, p.d)
# Wishart log pdf of R, using precomputed triangular L such that L*L' = R
Wishart_logpdf(R,L,p) = 0.5*(p.nu-p.d-1)*Lower.logdetsq(L,p.d) - 0.5*p.nu*matrix_dot(R,p.W,p.d) - p.c
Wishart_set_nu!(p,nu) = (Lower.scale!(p.M,sqrt(p.nu/nu),p.d); p.nu=nu; p.c=Wishart_constant(p))
Wishart_notify_W!(p) = (Lower.Cholesky_inverse!(p.M,p.A,p.W,p.d); Lower.scale!(p.M,sqrt(1/p.nu),p.d);
    p.c=Wishart_constant(p))
Wishart_clear!(p) = (fill!(p.sum_R,0.); p.sum_logdetR=0.; p.n=0)
Wishart_adjoin!(p,R,logdetR) = (for i=1:p.d*p.d; p.sum_R[i]+=R[i]; end; p.sum_logdetR+=logdetR; p.n+=1)
# Wishart_remove!(p) = ...


# # LogGamma distribution (X = log(Y) where Y~Gamma(a,b))
type LogGamma_params
    a::Float64
    b::Float64
    c::Float64
    LogGamma_params(a,b) = (p=new(); p.a=a; p.b=b; p.c=lgamma(a)-a*log(b); p)
end
LogGamma_logpdf(x,p) = p.a*x - p.b*exp(x) - p.c


type Hyperparameters
    m::MVN_params # parameters of m's (component means)
    R::Wishart_params # parameters of R's (component precisions)
    mm::MVN_params # parameters of m.m
    mR::Wishart_params # parameters of m.R
    Rn::LogGamma_params # parameters of R.nu
    RW::Wishart_params # parameters of R.W
    mp::MVN_params # temporary variable for the posterior of m's
    Rp::Wishart_params # temporary variable for the posterior of R's
    A::Array{Float64,2} # temporary variable for computations
    x::Array{Float64,1} # temporary variable for computations
    d::Int64 # dimension
end

# likelihood: Normal(x|mean=m,Cov=inv(R)) (Note: R is represented as L*L'.)
typealias Theta MVN_params
Theta_clear!,Theta_adjoin!,Theta_remove! = MVN_clear!,MVN_adjoin!,MVN_remove!
log_likelihood(x,p) = MVN_logpdf(x,p)
# prior: Normal(m|mean=H.m.m,Cov=inv(H.m.L*H.m.L')) Wishart(R|Scale=H.R.M*H.R.M',DOF=H.R.nu)
log_prior(p,H) = MVN_logpdf(p.m,H.m) + Wishart_logpdf(MVN_get_R!(p),p.L,H.R)
new_theta(H) = MVN_params(zeros(H.d),eye(H.d))
prior_sample!(p,H) = (MVN_sample!(p.m,H.m); Wishart_sample!(p.L,H.R); MVN_notify_L!(p))


sum_outer(x,m,d,n) = (S=zeros(d,d); for k = 1:n, i=1:d, j=1:d; S[i,j] += (x[k][i]-m[i])*(x[k][j]-m[j]); end; S)

function construct_hyperparameters(options)
    x = options.x
    # x should be an array of n data points of dimension d, e.g., an array of vectors, with n>=d
    n,d = length(x),length(x[1])
    m_hat = mean(x)
    C_hat = sum_outer(x,m_hat,d,n)/n
    R_hat = inv(cholfact(C_hat))
    
    mm = MVN_params(m_hat,R_hat) # these are fixed
    mR = Wishart_params(d,C_hat) # these are fixed
    
    m = MVN_params(m_hat,R_hat) # initialize to typical values (E(m.m), E(m.R))
    
    RW = Wishart_params(d,R_hat) # these are fixed
    
    # We take log(R.nu-d+1) ~ LogGamma(a,b). (Equivalently, R.nu-d+1 ~ Gamma(a,b).)
    # (This is different from Rasmussen, who takes 1/(R.nu-d+1) ~ Gamma(a,b).)
    Rn = LogGamma_params(2.,2.) # these are fixed
    
    R = Wishart_params(d,C_hat) # initialize to typical values
    # This choice of initial R.nu,R.W is derived by setting R.nu=E(R.nu)=d and R.W=E(R.W)=C_hat.
    
    mp = MVN_params(zeros(d),eye(d))
    Rp = Wishart_params(d,eye(d))
    
    A = zeros(d,d)
    x0 = zeros(d)
    
    return Hyperparameters(m,R,mm,mR,Rn,RW,mp,Rp,A,x0,d)
end


function NormalWishart_update!(a,b,m_params,R_params,H,active,density)
    d = H.d
    a_R = MVN_get_R!(a)
    m_params_R = MVN_get_R!(m_params)

    # update mean
    for i = 1:d*d; H.A[i] = m_params_R[i] + b.n * a_R[i]; end
    Lower.Cholesky!(H.mp.L, H.A, d); MVN_notify_L!(H.mp)
    for i = 1:d
        r = 0.; for j = 1:d; r += m_params_R[i,j]*m_params.m[j] + a_R[i,j]*b.sum_x[j]; end
        H.x[i] = r
    end
    Lower.solve_Lx_eq_y!(H.mp.L, H.x, H.x, d)
    Lower.solve_Ltx_eq_y!(H.mp.L, H.mp.m, H.x, d)
    if active; MVN_sample!(b.m, H.mp); end
    
    # update precision matrix
    Wishart_set_nu!(H.Rp, b.n + R_params.nu)
    for i = 1:d, j = 1:d
        H.Rp.W[i,j] = (R_params.nu*R_params.W[i,j] + b.sum_xx[i,j] 
            - b.sum_x[i]*b.m[j] - b.m[i]*b.sum_x[j] + b.n*b.m[i]*b.m[j]) / H.Rp.nu
    end; Wishart_notify_W!(H.Rp)
    if active; Wishart_sample!(b.L, H.Rp); MVN_notify_L!(b); end
    
    # compute density
    return (density? MVN_logpdf(b.m, H.mp) + Wishart_logpdf(MVN_get_R!(b), b.L, H.Rp) : NaN)
end

function Wishart_update!(p,W_params,nu_params,H)
    d = H.d
    
    # update W (inverse mean matrix)
    Wishart_set_nu!(H.Rp, W_params.nu + p.n * p.nu)
    for i = 1:d*d
        H.Rp.W[i] = (W_params.nu*W_params.W[i] + p.nu*p.sum_R[i]) / H.Rp.nu
    end; Wishart_notify_W!(H.Rp)
    Wishart_sample!(H.A, H.Rp) # H.A*H.A' ~ Wishart(Scale=H.Rp.M*H.Rp.M', DOF=H.Rp.nu)
    Lower.multiplyMNt!(p.W, H.A, H.A, d); Wishart_notify_W!(p) # p.W = H.A*H.A'
    
    # update nu (degrees of freedom) using a Metropolis move
    nu_0 = p.nu
    x_0 = log(p.nu-d+1)
    log_prior_0 = LogGamma_logpdf(x_0,nu_params)
    log_lik_0 = 0.5*(p.nu-d-1)*p.sum_logdetR - 0.5*p.nu*matrix_dot(p.W,p.sum_R,d) - p.n*p.c
    
    # make proposal
    x_prop = x_0 + randn()*NU_SIGMA_PROP
    Wishart_set_nu!(p, exp(x_prop)+d-1) # Note: p.c is updated here
    log_prior_prop = LogGamma_logpdf(x_prop,nu_params)
    log_lik_prop = 0.5*(p.nu-d-1)*p.sum_logdetR - 0.5*p.nu*matrix_dot(p.W,p.sum_R,d) - p.n*p.c
    
    # accept or reject
    p_accept = min(1.0, exp(log_prior_prop + log_lik_prop - log_prior_0 - log_lik_0))
    if rand()>p_accept
        Wishart_set_nu!(p,nu_0) # if we don't accept, set nu back to it's original value
    end
end


function update_parameter!(theta_a,theta_b,H,active,density)
    return NormalWishart_update!(theta_a,theta_b,H.m,H.R,H,active,density)
end

function update_hyperparameters!(H,theta,list,t,x,z)
    # update the hyperparameters of the mean (m.m and m.R)
    MVN_clear!(H.m)
    for k = 1:t; MVN_adjoin!(H.m, theta[list[k]].m); end
    NormalWishart_update!(H.m,H.m,H.mm,H.mR,H,true,false)
    
    # update the hyperparameters of the precision (R.W and R.nu)
    Wishart_clear!(H.R)
    for k = 1:t; tk=theta[list[k]]; Wishart_adjoin!(H.R, MVN_get_R!(tk), tk.logdetR); end
    Wishart_update!(H.R,H.RW,H.Rn,H)
end



# sample from a discrete distribution with weights proportional to q = [q_1,...,q_n]
randp(q,n) = (cdf=cumsum(q/sum(q)); [findfirst(rand().<cdf)::Int64 for i = 1:n])

# sample from a mixture
function mixrnd(n,probabilities,theta)
    k = length(probabilities)
    z = randp(probabilities[1:k],n) # component assignments
    d = theta[1].d
    x = [MVN_sample!(zeros(d),theta[z[i]])::Array{Float64,1} for i = 1:n]
    return x,z
end

function mixture_density(x,probabilities,theta)
    k = length(probabilities)
    px = 0.0
    for i = 1:k
        px += probabilities[i]*exp(log_likelihood(x,theta[i]))
    end
    return px
end

end # module MVNmodel
using .MVNmodel

# Include generic code
include("generic.jl")

# Include core sampler code
include("coreNonconjugate.jl")

end # module MVN







