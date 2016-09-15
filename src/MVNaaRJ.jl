# Reversible jump MCMC for axis-aligned multivariate normal mixtures, using a conjugate prior.
module MVNaaRJ

using Distributions

typealias Data Array{Float64,1}
type Theta
    mu::Array{Float64,1}     # means
    lambda::Array{Float64,1} # precisions
end

# Include generic code
include("generic.jl")

logsumexp(a,b) = (m = max(a,b); m == -Inf ? -Inf : log(exp(a-m) + exp(b-m)) + m)

lognormpdf(x,m,v) = -0.5*(x-m)*(x-m)/v - 0.5*log(2*pi*v)
logmvnpdf(x,m,v,c) = (l=0.0; for j=1:length(x); l += lognormpdf(x[j],m[c,j],v[c,j]); end; l)

type Hyperparameters
    log_pk_vec::Array{Float64,1} # log of prior on the number of components
    gamma::Float64 # Dirichlet parameter
    d::Int64    # dimension
    # m[c,i]|v[c,i] ~ N(m,v[c,i]/r)
    m::Float64  # prior mean of m's
    r::Float64  # factor for prior precision of m's
    # v[c,i] ~ InvGamma(a,b)
    a::Float64  # prior shape of v's
    b::Float64  # prior scale of v's
end

function construct_hyperparameters(options)
    x = options.x
    n = length(x)
    d = length(x[1])
    mu = mean(x)
    v = mean([xi.*xi for xi in x]) - mu.*mu  # sample variance
    @assert(all(abs(mu) .< 1e-10) && all(abs(v - 1.0) .< 1e-10), "Data must be normalized to zero mean, unit variance.")
    m = 0.0
    r = 1.0
    a = 1.0
    b = 1.0
    gamma = options.gamma
    log_pk_fn = eval(parse(options.log_pk))
    log_pk_vec = Float64[log_pk_fn(k) for k = 1:options.k_max]
    return Hyperparameters(log_pk_vec,gamma,d,m,r,a,b)
end


type List
    first::Int64
    next::Array{Int64,1}
    List(maxlen) = (l=new(); l.first=0; l.next=zeros(Int64,maxlen); return l)
end
insert!(l::List,i) = (l.next[i] = l.first; l.first = i) # Note: i must not already be in the list.
pop!(l::List) = (i = l.first; @assert(i!=0,"list capacity exceeded"); l.first = l.next[i]; i)
function remove!(l::List,i) # Note: i must already be in the list.
    if i==l.first
        l.first = l.next[i]
    else
        j = l.first
        while j != 0
            if i==l.next[j]; l.next[j] = l.next[i]; break; end
            j = l.next[j]
        end
    end
end
#pr(l::List) = (i=l.first; while i!=0; print(i," "); i=l.next[i]; end; println())

type IndexList
    first::Array{Int64,1}
    next::Array{Int64,1}
    IndexList(kmax,n) = (l=new(); l.first=zeros(Int64,kmax); l.next=zeros(Int64,n); return l)
end
insert!(l::IndexList,c,i) = (l.next[i] = l.first[c]; l.first[c] = i)
clear!(l::IndexList,c) = (l.first[c] = 0)


# sample from categorical distn p on list
function randp(p,list::List)
    u = rand()
    i = list.first
    s = p[i]
    while u > s
        i = list.next[i]
        s += p[i]
    end
    return i
end

# sample from uniform distn on list of k items
function randunif(k,list::List)
    u = rand()*k
    i = list.first
    s = 1.0
    while u > s
        i = list.next[i]
        s += 1.0
    end
    return i
end

# update assignments z
function update_assignments!(x,z,w,m,v,clist,log_p,p)
    for i = 1:length(x)
        c,log_s = clist.first,-Inf
        while c != 0
            log_p[c] = log(w[c]) + logmvnpdf(x[i],m,v,c)
            log_s = logsumexp(log_s,log_p[c])
            c = clist.next[c]
        end
        c = clist.first
        while c != 0
            p[c] = exp(log_p[c] - log_s)
            c = clist.next[c]
        end
        z[i] = randp(p,clist)
    end
end

# recompute t,s,sx,ssx,ilist
function recompute_statistics!(x,z,s,sx,ssx,ilist,clist,H)
    c = clist.first
    while c != 0
        s[c] = 0
        clear!(ilist,c)
        for j = 1:H.d; sx[c,j] = 0.0; ssx[c,j] = 0.0; end
        c = clist.next[c]
    end
    t = 0
    for i = 1:length(x)
        c = z[i]
        if s[c]==0; t+=1; end
        insert!(ilist,c,i)
        s[c] += 1
        for j = 1:H.d
            sx[c,j] += x[i][j]
            ssx[c,j] += x[i][j]*x[i][j]
        end
    end
    return t
end

# update weights, means, and variances
function update_parameters!(w,m,v,s,sx,ssx,clist,H)
    # update weights (unnormalized for now)
    c = clist.first
    w_sum = 0.0
    while c != 0
        w[c] = rand(Gamma(H.gamma+s[c],1))
        w_sum += w[c]
        c = clist.next[c]
    end
    
    c = clist.first
    while c != 0
        A = H.a + 0.5*(s[c]+1)
        for j = 1:H.d
            # update mean
            L = (H.r + s[c])/v[c,j]
            M = (H.r*H.m + sx[c,j]) / (H.r + s[c])
            m[c,j] = randn()/sqrt(L) + M
            
            # update variance
            B = H.b + 0.5*(ssx[c,j] - 2*m[c,j]*sx[c,j] + s[c]*m[c,j]*m[c,j]) + 0.5*H.r*(m[c,j]-H.m)*(m[c,j]-H.m)
            # if (A<=0 || B<=0); println(A," ",B)
                # println(H.b)
                # println(0.5*(ssx[c,j] - 2*m[c,j]*sx[c,j] + s[c]*m[c,j]*m[c,j]))
                # println(0.5*H.r*(m[c,j]-H.m)*(m[c,j]-H.m))
            # end
            v[c,j] = rand(InverseGamma(A,B))
        end
        
        # normalize weights
        w[c] = w[c]/w_sum
        c = clist.next[c]
    end
end

function log_prior_ratio(k,w,m,v,cm,c1,c2,s1,s2,H)
    ssm,slv,siv = 0.0,0.0,0.0
    for j = 1:H.d
        ssm += (m[c1,j]-H.m)^2/v[c1,j] + (m[c2,j]-H.m)^2/v[c2,j] - (m[cm,j]-H.m)^2/v[cm,j]
        slv += log(v[c1,j]*v[c2,j]/v[cm,j])
        siv += 1/v[c1,j] + 1/v[c2,j] - 1/v[cm,j]
    end
    log_rprior = (H.log_pk_vec[k+1] - H.log_pk_vec[k] # prior on k
               + (H.gamma-1)*log(w[c1]*w[c2]/w[cm]) - lbeta(H.gamma,k*H.gamma) # prior on w
               + s1*log(w[c1]) + s2*log(w[c2]) - (s1+s2)*log(w[cm]) # prior on z
               + 0.5*H.d*log(H.r/(2*pi)) - 0.5*H.r*ssm - 0.5*slv # prior on m
               + H.d*H.a*log(H.b) - H.d*lgamma(H.a) + (-H.a-1)*slv - H.b*siv) # prior on v
    return log_rprior
end

function log_Jacobian(w,v,u1,u2,cm,H)
    log_J = log(w[cm]) - 1.5*H.d*log(u1*(1-u1))
    for j = 1:H.d; log_J += log(1-u2[j]*u2[j]) + 1.5*log(v[cm,j]); end
    return log_J
end

# propose split-merge move
function split_merge!(k,t,x,z,w,m,v,s,sx,ssx,clist,cfree,ilist, zs,u2,u3,du1,du2,du3,H)
    if ((k==1) || (rand()<0.5)) # propose split
        # choose component indices
        cm = randunif(k,clist)
        c1 = pop!(cfree)
        c2 = pop!(cfree)
        # sample u's
        log_pu = 0.0
        u1 = rand(du1); log_pu += logpdf(du1,u1)
        for j = 1:H.d
            u2[j] = rand(du2); log_pu += logpdf(du2,u2[j])
            u3[j] = rand(du3); log_pu += logpdf(du3,u3[j])
        end
        # split w, m, and v
        w[c1],w[c2] = w[cm]*u1, w[cm]*(1-u1)
        for j = 1:H.d
            m[c1,j] = m[cm,j] - u2[j]*sqrt(v[cm,j]*w[c2]/w[c1])
            m[c2,j] = m[cm,j] + u2[j]*sqrt(v[cm,j]*w[c1]/w[c2])
            v[c1,j] = u3[j]*(1-u2[j]*u2[j])*v[cm,j]*w[cm]/w[c1]
            v[c2,j] = (1-u3[j])*(1-u2[j]*u2[j])*v[cm,j]*w[cm]/w[c2]
        end
        # split z and simultaneously compute the log likelihood ratio
        log_pzs = 0.0
        log_rlik = 0.0
        s1,s2 = 0,0
        i = ilist.first[cm]
        while i != 0
            l1,l2,lm = logmvnpdf(x[i],m,v,c1),logmvnpdf(x[i],m,v,c2),logmvnpdf(x[i],m,v,cm)
            p1 = 1/(1 + exp(l2-l1)*w[c2]/w[c1])
            if rand() < p1
                zs[i] = c1
                s1 += 1
                log_pzs += log(p1)
                log_rlik += l1 - lm
            else
                zs[i] = c2
                s2 += 1
                log_pzs += log(1-p1)
                log_rlik += l2 - lm
            end
            i = ilist.next[i]
        end
        # log of the ratio of move probabilities (choice of split or merge)
        log_rmove = (k==1 ? log(0.5) : 0.0)
        
        # log prior ratio
        log_rprior = log_prior_ratio(k,w,m,v,cm,c1,c2,s1,s2,H)

        # log of the Jacobian
        log_J = log_Jacobian(w,v,u1,u2,cm,H)

        # acceptance probability
        p_accept = min(1,exp(log_rlik + log_rprior + log_rmove - log_pu - log_pzs + log_J))

        # accept or reject
        if rand() < p_accept
            # k,t,z,s,w,m,v,sx,ssx,clist,cfree,ilist
            k += 1
            t += Int((s1>0) && (s2>0))
            s[c1],s[c2] = s1,s2
            for j=1:H.d; sx[c1,j]=0.0; sx[c2,j]=0.0; ssx[c1,j]=0.0; ssx[c2,j]=0.0; end
            clear!(ilist,c1); clear!(ilist,c2)
            i = ilist.first[cm]
            while i != 0
                c = zs[i]
                z[i] = c
                for j = 1:H.d
                    sx[c,j] += x[i][j]
                    ssx[c,j] += x[i][j]*x[i][j]
                end
                next = ilist.next[i]
                insert!(ilist,c,i)
                i = next
            end
            insert!(clist,c1); insert!(clist,c2)
            remove!(clist,cm); insert!(cfree,cm)
        else
            insert!(cfree,c1); insert!(cfree,c2)
        end
    else # propose merge
        # choose component indices
        c1 = randunif(k,clist); remove!(clist,c1)
        c2 = randunif(k-1,clist); insert!(clist,c1)
        cm = pop!(cfree)
        # merge w, m, and v
        w[cm] = w[c1] + w[c2]
        for j = 1:H.d
            m[cm,j] = (w[c1]*m[c1,j] + w[c2]*m[c2,j])/w[cm]
            v[cm,j] = (w[c1]*(m[c1,j]^2 + v[c1,j]) + w[c2]*(m[c2,j]^2 + v[c2,j]))/w[cm] - m[cm,j]^2
        end
        # compute u's
        log_pu = 0.0
        u1 = w[c1]/w[cm]; log_pu += logpdf(du1,u1)
        for j = 1:H.d
            u2[j] = (m[cm,j] - m[c1,j])/sqrt(v[cm,j]*w[c2]/w[c1]); log_pu += logpdf(du2,u2[j])
            u3[j] = (v[c1,j]*w[c1])/(v[cm,j]*w[cm]*(1-u2[j]*u2[j])); log_pu += logpdf(du3,u3[j])
        end
        # merge z and compute the log likelihood ratio
        log_pzs = 0.0
        log_rlik = 0.0
        i = ilist.first[c1]
        while i != 0
            l1,l2,lm = logmvnpdf(x[i],m,v,c1),logmvnpdf(x[i],m,v,c2),logmvnpdf(x[i],m,v,cm)
            p1 = 1/(1 + exp(l2-l1)*w[c2]/w[c1])
            log_pzs += log(p1)
            log_rlik += lm - l1
            i = ilist.next[i]
        end
        i = ilist.first[c2]
        while i != 0
            l1,l2,lm = logmvnpdf(x[i],m,v,c1),logmvnpdf(x[i],m,v,c2),logmvnpdf(x[i],m,v,cm)
            p1 = 1/(1 + exp(l2-l1)*w[c2]/w[c1])
            log_pzs += log(1-p1)
            log_rlik += lm - l2
            i = ilist.next[i]
        end
        # log of the ratio of move probabilities (choice of split or merge)
        log_rmove = 0.0

        # log prior ratio
        log_rprior = -log_prior_ratio(k-1,w,m,v,cm,c1,c2,s[c1],s[c2],H)

        # log of the Jacobian
        log_J = -log_Jacobian(w,v,u1,u2,cm,H)
        
        # acceptance probability
        p_accept = min(1,exp(log_rlik + log_rprior + log_rmove + log_pu + log_pzs + log_J))

        # accept or reject
        if rand() < p_accept
            # k,t,z,w,m,v,s,sx,ssx,clist,cfree,ilist
            k -= 1
            t -= Int((s[c1]>0) && (s[c2]>0))
            s[cm] = s[c1] + s[c2]
            for j = 1:H.d
                sx[cm,j] = sx[c1,j] + sx[c2,j]
                ssx[cm,j] = ssx[c1,j] + ssx[c2,j]
            end
            clear!(ilist,cm)
            for i in [ilist.first[c1],ilist.first[c2]]
                while i != 0
                    z[i] = cm
                    next = ilist.next[i]
                    insert!(ilist,cm,i)
                    i = next
                end
            end
            insert!(clist,cm) 
            remove!(clist,c1); remove!(clist,c2)
            insert!(cfree,c1); insert!(cfree,c2)
        else
            insert!(cfree,cm)
        end
    end
    return k,t
end

function birth_death!(k,t,n,w,m,v,s,sx,ssx,clist,cfree,ilist,H)
    k0 = k-t # number of empty components
    dv = InverseGamma(H.a,H.b)
    if rand() < 0.5 # propose birth of an empty component
        c = pop!(cfree)
        w[c] = rand(Beta(1,k))
        for j = 1:H.d
            v[c,j] = rand(dv)
            m[c,j] = randn()*sqrt(v[c,j]/H.r) + H.m
        end
        log_A = (H.log_pk_vec[k+1] - H.log_pk_vec[k]
              + (n+k*(H.gamma-1))*log(1-w[c]) + (H.gamma-1)*log(w[c])
              + lbeta(1,k) - lbeta(H.gamma,k*H.gamma) + log(k+1) - log(k0+1))
        p_accept = min(1,exp(log_A))
        
        # accept or reject
        if rand() < p_accept
            # k,t,z,w,m,v,s,sx,ssx,clist,cfree,ilist
            k += 1
            s[c] = 0
            for j = 1:H.d; sx[c,j] = 0.0; ssx[c,j] = 0.0; end
            insert!(clist,c)
            clear!(ilist,c)
        else
            insert!(cfree,c)
        end
    elseif (k0 > 0) # propose death of an empty component
        # randomly choose an empty component
        u = rand()*k0
        c = clist.first
        S = Int(s[c]==0)
        while u > S
            c = clist.next[c]
            S += Int(s[c]==0)
        end
        # compute acceptance probability
        log_A = (H.log_pk_vec[k] - H.log_pk_vec[k-1]
              + (n+(k-1)*(H.gamma-1))*log(1-w[c]) + (H.gamma-1)*log(w[c])
              + lbeta(1,k-1) - lbeta(H.gamma,(k-1)*H.gamma) + log(k) - log(k0))
        p_accept = min(1,exp(-log_A))
        
        # accept or reject
        if rand() < p_accept
            # k,t,z,w,m,v,s,sx,ssx,clist,cfree,ilist
            k -= 1
            remove!(clist,c)
            insert!(cfree,c)
        end
    end
    return k,t
end


function sampler(options,n_total,n_keep)
    H = construct_hyperparameters(options)
    x,n = options.x,options.n
    kmax = min(n,options.k_max)
    @assert(n==length(x))
    d = H.d # dimension
    du1,du2,du3 = Beta(2,2),Beta(2,2),Beta(1,1) # distns of aux vars u1,u2,u3
    
    # State
    # k,t,z,s,w,m,v,sx,ssx,cfirst,cnext,cfree,ilist
    k = 1 # number of components
    t = 1 # number of clusters (non-empty components)
    z = zeros(Int64,n) # z[i] = component assignment of datapoint i
    w = zeros(kmax) # w[c] = mixture weight of component c
    m = zeros(kmax,d) # m[c,:] = means of component c
    v = zeros(kmax,d) # v[c,:] = variances of component c
    s = zeros(Int64,kmax) # s[c] = # of points assigned to component c
    sx = zeros(kmax,d) # sx[c,:] = sum of the x's in component c
    ssx = zeros(kmax,d) # ssx[c,:] = sum of the squares of the x's in component c
    clist = List(kmax) # component list
    cfree = List(kmax) # list of unused components
    ilist = IndexList(kmax,n) # lists of indices assigned to each component
    
    # initialize with a single component
    z[:] = 1
    s[1] = n
    m[1,:] = mean(x) # initialize to the sample mean
    v[1,:] = mean([(xi-vec(m[1,:])).^2 for xi in x]) # initialize to the sample variance
    sx[1,:] = sum(x)
    ssx[1,:] = sum([xi.^2 for xi in x])
    insert!(clist,1)
    for c=kmax:-1:2; insert!(cfree,c); end
    for i = 1:n; insert!(ilist,1,i); end
    
    # temporary variables
    p = zeros(kmax) # temporary variable for randp
    log_p = zeros(kmax) # temporary variable
    zs = zeros(Int64,n) # temporary assignments for split-merge
    u2,u3 = zeros(d),zeros(d) # temporary variables for split-merge

    # iterations to record
    keepers = round(Int,linspace(round(Int,n_total/n_keep),n_total,n_keep))
    keep_index = 0
    
    # records
    #k_r = zeros(Int64,n_total)
    t_r = zeros(Int8,n_total); @assert(kmax < 2^7)
    N_r = zeros(Int16,kmax+3,n_total); @assert(n < 2^15)
    z_r = zeros(Int8,n,n_keep); @assert(kmax < 2^7)
    theta_r = Array(Theta,kmax+3,n_keep)
    
    for iteration = 1:n_total
        # update weights, means, and variances
        update_parameters!(w,m,v,s,sx,ssx,clist,H)

        # update assignments z
        update_assignments!(x,z,w,m,v,clist,log_p,p)
        # recompute t,s,sx,ssx,ilist based on the new z
        t = recompute_statistics!(x,z,s,sx,ssx,ilist,clist,H)
        
        # split/merge
        k,t = split_merge!(k,t,x,z,w,m,v,s,sx,ssx,clist,cfree,ilist, zs,u2,u3,du1,du2,du3,H)

        # birth/death
        k,t = birth_death!(k,t,n,w,m,v,s,sx,ssx,clist,cfree,ilist,H)

        # record-keeping
        #k_r[iteration] = k
        t_r[iteration] = t
        for i = 1:n
            N_r[z[i],iteration] += 1
        end
        if iteration==keepers[keep_index+1]
            keep_index += 1
            for i = 1:n; z_r[i,keep_index] = z[i]; end
            c = clist.first
            while c != 0
                theta_r[c,keep_index] = Theta(m[c,:][:], 1.0./v[c,:][:])
                c = clist.next[c]
            end
        end
    end
    
    return t_r,N_r,z_r,theta_r,keepers
end

end # module



















