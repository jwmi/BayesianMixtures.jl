
logsumexp(a,b) = (m = max(a,b); m == -Inf ? -Inf : log(exp(a-m) + exp(b-m)) + m)

function randp(p,k)
    s = 0.; for j = 1:k; s += p[j]; end
    u = rand()*s
    j = 1
    C = p[1]
    while u > C
        j += 1
        C += p[j]
    end
    @assert(j <= k)
    return j
end

function randlogp!(log_p,k)
    log_s = -Inf; for j = 1:k; log_s = logsumexp(log_s,log_p[j]); end
    p = log_p
    for j = 1:k; p[j] = exp(log_p[j]-log_s); end
    return randp(p,k)
end

function ordered_insert!(index,list,t)
    j = t
    while (j>0) && (list[j]>index)
        list[j+1] = list[j]
        j -= 1
    end
    list[j+1] = index
end

function ordered_remove!(index,list,t)
    for j = 1:t
        if list[j]>=index; list[j] = list[j+1]; end
    end
end

function ordered_next(list)
    j = 1
    while list[j]==j; j += 1; end
    return j
end
                

function restricted_Gibbs!(zsa,zsb,tia,tib,tja,tjb,cia,cib,cja,cjb,ni,nj,i,j,S,ns,x,b,H,active)
# NOTE: The sufficient statistics of tia and tja must be in sync with zsa.
# Also, note that the sufficient statistics of tib and tjb are not updated in this procedure.
    log_p = 0.
    for ks = 1:ns
        k = S[ks]
        if k!=i && k!=j
            if zsa[k]==cia; ni -= 1; else; nj -= 1; end
            Li = log_likelihood(x[k],tia)
            Lj = log_likelihood(x[k],tja)
            Pi = exp(log(ni+b)+Li - logsumexp(log(ni+b)+Li,log(nj+b)+Lj))
            
            if active
                if rand()<Pi
                    if zsa[k]==cja; Theta_remove!(tja,x[k]); Theta_adjoin!(tia,x[k]); end
                    zsb[k] = cib
                else
                    if zsa[k]==cia; Theta_remove!(tia,x[k]); Theta_adjoin!(tja,x[k]); end
                    zsb[k] = cjb
                end
            end
            if zsb[k]==cib
                ni += 1
                log_p += log(Pi)
            else
                nj += 1
                log_p += log(1-Pi)
            end
        end
    end
    log_p += update_parameter!(tia,tib,H,active,true)
    log_p += update_parameter!(tja,tjb,H,active,true)
    return log_p,ni,nj
end

function split_merge!(x,z,zs,S,theta,list,N,t,H,a,b,log_v,n_split,n_merge)
    n = length(x)
    
    # randomly choose a pair of indices
    i = round(Int,ceil(rand()*n))
    j = round(Int,ceil(rand()*(n-1))); if j>=i; j += 1; end
    ci0,cj0 = z[i],z[j]
    ti0,tj0 = theta[ci0],theta[cj0]
    
    # set S[1],...,S[ns] to the indices of the points in clusters ci0 and cj0
    ns = 0
    for k = 1:n
        if z[k]==ci0 || z[k]==cj0; ns += 1; S[ns] = k; end
    end
    
    # find available cluster IDs for merge and split parameters
    k = 1
    while list[k]==k; k += 1; end; cm = k
    while list[k]==k+1; k += 1; end; ci = k+1
    while list[k]==k+2; k += 1; end; cj = k+2
    tm,ti,tj = theta[cm],theta[ci],theta[cj]
    
    # randomly choose the merge launch state
    prior_sample!(tm,H)  # sample initial parameter from the prior
    for ks = 1:ns; Theta_adjoin!(tm,x[S[ks]]); end # get the sufficient statistics
    for rep = 1:n_merge  # make several moves
        update_parameter!(tm,tm,H,true,false)
    end
    
    # randomly choose the split launch state
    prior_sample!(ti,H)  # sample initial parameters from the prior
    prior_sample!(tj,H)
    zs[i] = ci; Theta_adjoin!(ti,x[i]); ni = 1
    zs[j] = cj; Theta_adjoin!(tj,x[j]); nj = 1
    for ks = 1:ns  # start with a uniformly chosen split
        k = S[ks]
        if k!=i && k!=j
            if rand()<0.5; zs[k] = ci; Theta_adjoin!(ti,x[k]); ni += 1
            else;          zs[k] = cj; Theta_adjoin!(tj,x[k]); nj += 1
            end
        end
    end
    for rep = 1:n_split  # make several moves
        log_p,ni,nj = restricted_Gibbs!(zs,zs,ti,ti,tj,tj,ci,ci,cj,cj,ni,nj,i,j,S,ns,x,b,H,true)
    end
    
    # make proposal
    if ci0==cj0  # propose a split
        # make one final sweep and compute it's probability density
        log_prop_ab,ni,nj = restricted_Gibbs!(zs,zs,ti,ti,tj,tj,ci,ci,cj,cj,ni,nj,i,j,S,ns,x,b,H,true)
        
        # compute probability density of going from merge launch state to original state
        log_prop_ba = update_parameter!(tm,ti0,H,false,true)
        
        # compute acceptance probability
        log_prior_b = log_v[t+1] + lgamma(ni+b)+lgamma(nj+b)-2*lgamma(a) + log_prior(ti,H) + log_prior(tj,H)
        log_prior_a = log_v[t] + lgamma(ns+b)-lgamma(a) + log_prior(ti0,H)
        log_lik_ratio = 0.
        for ks = 1:ns; k = S[ks]
            log_lik_ratio += log_likelihood(x[k],(zs[k]==ci? ti:tj)) - log_likelihood(x[k],ti0)
        end
        p_accept = min(1.0, exp(log_prop_ba-log_prop_ab + log_prior_b-log_prior_a + log_lik_ratio))
        #println("split proposal: ",p_accept)
        
        # accept or reject
        if rand()<p_accept # accept split
            #  z, list, N, theta, t
            for ks = 1:ns; z[S[ks]] = zs[S[ks]]; end
            ordered_remove!(ci0,list,t)
            ordered_insert!(ci,list,t-1)
            ordered_insert!(cj,list,t)
            N[ci0],N[ci],N[cj] = 0,ni,nj
            t += 1
            Theta_clear!(ti0)
        else # reject split
            Theta_clear!(ti)
            Theta_clear!(tj)
        end
        Theta_clear!(tm)
        
    else  # propose a merge
        # make one final sweep and compute its probability density
        log_prop_ab = update_parameter!(tm,tm,H,true,true)
        
        # compute probability density of going from split launch state to original state
        log_prop_ba,ni,nj = restricted_Gibbs!(zs,z,ti,ti0,tj,tj0,ci,ci0,cj,cj0,ni,nj,i,j,S,ns,x,b,H,false)
        
        # compute acceptance probability
        log_prior_b = log_v[t-1] + lgamma(ns+b)-lgamma(a) + log_prior(tm,H)
        log_prior_a = log_v[t] + lgamma(ni+b)+lgamma(nj+b)-2*lgamma(a) + log_prior(ti0,H) + log_prior(tj0,H)
        log_lik_ratio = 0.
        for ks = 1:ns; k = S[ks]
            log_lik_ratio += log_likelihood(x[k],tm) - log_likelihood(x[k],(z[k]==ci0? ti0:tj0))
        end
        p_accept = min(1.0, exp(log_prop_ba-log_prop_ab + log_prior_b-log_prior_a + log_lik_ratio))
        #println("merge proposal: ",p_accept)
        
        # accept or reject
        if rand()<p_accept # accept merge
            #  z, list, N, theta, t
            for ks = 1:ns; z[S[ks]] = cm; end
            ordered_remove!(ci0,list,t)
            ordered_remove!(cj0,list,t-1)
            ordered_insert!(cm,list,t-2)
            N[cm],N[ci0],N[cj0] = ns,0,0
            t -= 1
            Theta_clear!(ti0)
            Theta_clear!(tj0)
        else # reject merge
            Theta_clear!(tm)
        end
        Theta_clear!(ti)
        Theta_clear!(tj)
    end
    return t
end

function sampler(options,n_total,n_keep)
    x,n = options.x,options.n
    t_max = options.t_max
    a,b,log_v = options.a,options.b,options.log_v
    use_splitmerge,n_split,n_merge = options.use_splitmerge,options.n_split,options.n_merge
    use_hyperprior = options.use_hyperprior

    model_type = options.model_type
    alpha_random,alpha = options.alpha_random,options.alpha
    sigma_alpha = 0.1 # scale for MH proposals in alpha move

    @assert(n==length(x))
    keepers = zeros(Int,n_keep)
    keepers[:] = round.(Int,linspace(round(Int,n_total/n_keep),n_total,n_keep))
    keep_index = 0

    t = 1  # number of clusters
    z = ones(Int,n)  # z[i] = the cluster ID associated with data point i
    list = zeros(Int,t_max+3); list[1] = 1  # list[1:t] = the list of active cluster IDs
                                            # list is maintained in increasing order for 1:t, and is 0 after that.
                                            # Further, maximum(list[1:t]) is always <= t_max+3
    c_next = 2  # an available cluster ID to be used next
    N = zeros(Int,t_max+3); N[1] = n  # N[c] = size of cluster c

    H = construct_hyperparameters(options)
    theta = [new_theta(H)::Theta for c = 1:t_max+3]  # theta[c] = parameters for cluster c

    prior_sample!(theta[1],H)  # sample parameters for the initial cluster
    for i = 1:n; Theta_adjoin!(theta[1],x[i]); end

    log_p = zeros(n+1)
    zs = ones(Int,n)  # temporary variable used for split-merge assignments
    S = zeros(Int,n)  # temporary variable used for split-merge indices
    
    log_Nb = log.((1:n) + b)
    
    # Record-keeping variables
    t_r = zeros(Int8,n_total); @assert(t_max < 2^7)
    N_r = zeros(Int16,t_max+3,n_total); @assert(n < 2^15)
    z_r = zeros(Int8,n,n_keep); @assert(t_max < 2^7)
    theta_r = Array{Theta}(t_max+3,n_keep)
    
    for iteration = 1:n_total
        #  -------------- Resample thetas and H --------------
        for j = 1:t; c = list[j]
            update_parameter!(theta[c],theta[c],H,true,false)
        end
        if use_hyperprior
            update_hyperparameters!(H,theta,list,t,x,z)
        end
        if model_type=="DPM" && alpha_random
            # Metropolis-Hastings move for DP concentration parameter (using p_alpha(a) = exp(-a) = Exp(a|1))
            aprop = alpha*exp(randn()*sigma_alpha)
            top = t*log(aprop) - lgamma(aprop+n) + lgamma(aprop) - aprop + log(aprop)
            bot = t*log(alpha) - lgamma(alpha+n) + lgamma(alpha) - alpha + log(alpha)
            if rand() < min(1.0,exp(top-bot))
                alpha = aprop
            end
            #for i=1:t_max+1; log_v[i] = i*log(alpha) - lgamma(alpha+n) + lgamma(alpha); end
            log_v = float(1:t_max+1)*log(alpha) - lgamma(alpha+n) + lgamma(alpha)
        end
        
        # -------------- Resample z's --------------
        for i = 1:n
            # remove point i from it's cluster
            c = z[i]    
            N[c] -= 1
            if N[c]>0
                c_prop = c_next
                prior_sample!(theta[c_prop],H)
            else
                c_prop = c
                # remove cluster {i}, keeping list in proper order
                ordered_remove!(c,list,t)
                t -= 1
            end
            
            # compute probabilities for resampling
            for j = 1:t; cc = list[j]
                log_p[j] = log_Nb[N[cc]] + log_likelihood(x[i],theta[cc])
            end
            log_p[t+1] = log_v[t+1]-log_v[t] + log(a) + log_likelihood(x[i],theta[c_prop])
            
            # sample a new cluster for it
            j = randlogp!(log_p,t+1)
            
            # add point i to it's new cluster
            if j<=t
                c = list[j]
            else
                c = c_prop
                ordered_insert!(c,list,t)
                t += 1
                c_next = ordered_next(list)
                @assert(t<=t_max, "Sampled t has exceeded t_max. Increase t_max and retry.")
            end
            # update sufficient statistics if point i has changed clusters
            if c != z[i]
                Theta_remove!(theta[z[i]],x[i])
                Theta_adjoin!(theta[c],x[i])
            end
            z[i] = c
            N[c] += 1
        end
        
        
        # -------------- Split/merge move --------------
        if use_splitmerge
            t = split_merge!(x,z,zs,S,theta,list,N,t,H,a,b,log_v,n_split,n_merge)
            c_next = ordered_next(list)
            @assert(t<=t_max, "Sampled t has exceeded t_max. Increase t_max and retry.")
        end 
    
    
        # -------------- Record results --------------
        t_r[iteration] = t
        for j = 1:t
            N_r[list[j],iteration] = N[list[j]]
        end
        if iteration==keepers[keep_index+1]
            keep_index += 1
            for i = 1:n; z_r[i,keep_index] = z[i]; end
            for j = 1:t
                c = list[j]
                theta_r[c,keep_index] = deepcopy(theta[c])
            end
        end
    end
    
    return t_r,N_r,z_r,theta_r,keepers
end

















