
function randp(p,d)
    s = 0.; for j = 1:d; s += p[j]; end
    u = rand()*s
    j = 1
    C = p[1]
    while u > C
        j += 1
        C += p[j]
    end
    @assert(j <= d)
    return j
end

function ordered_insert_next!(list,t)
    j = t
    while (j>0) && (list[j]>j)
        list[j+1] = list[j]
        j -= 1
    end
    list[j+1] = j+1
    return j+1
end

function ordered_remove!(index,list,t)
    for j = 1:t
        if list[j]>=index; list[j] = list[j+1]; end
    end
end

function restricted_Gibbs!(zsa,zsb,ni,nj,xs,ns,is,js,tia,tib,tja,tjb,b,H,active)
    log_p = 0.
    for k = 1:ns
        if k!=is && k!=js
            if zsa[k]; ni -= 1; else; nj -= 1; end
            Li = likelihood(xs[k],tia)
            Lj = likelihood(xs[k],tja)
            ps = (ni+b)*Li/((ni+b)*Li + (nj+b)*Lj)
            if active; zsb[k] = (rand()<ps); end
            if zsb[k]
                ni += 1
                log_p += log(ps)
            else
                nj += 1
                log_p += log(1-ps)
            end
        end
    end
    
    log_p += update_parameter!(tia,tib,xs,zsb,true,H,active)
    log_p += update_parameter!(tja,tjb,xs,zsb,false,H,active)
    return log_p,ni,nj
end


function split_merge(xs,zs0,is,js,ti,tj,tm,ti0,tj0,t,H,a,b,log_v,n_split,n_merge)
    ns = length(xs)
    splitting = zs0[js]
    dummy = zeros(Int,ns)
    
    # randomly choose the split launch state
    zs = (rand(ns).<0.5)  # start with a uniformly chosen split
    zs[is],zs[js] = true,false
    ni = sum(zs); nj = ns-ni  # number of points in the two clusters
    prior_sample!(ti,H)  # sample initial parameters from the prior
    prior_sample!(tj,H)
    for rep = 1:n_split  # make several moves
        log_p,ni,nj = restricted_Gibbs!(zs,zs,ni,nj,xs,ns,is,js,ti,ti,tj,tj,b,H,true)
    end
    
    # randomly choose the merge launch state
    prior_sample!(tm,H)  # sample initial parameter from the prior
    for rep = 1:n_merge  # make several moves
        update_parameter!(tm,tm,xs,dummy,0,H,true)
    end
    
    # make proposal
    if splitting  # propose a split
        # make one final sweep and compute it's probability density
        log_prop_ab,ni,nj = restricted_Gibbs!(zs,zs,ni,nj,xs,ns,is,js,ti,ti,tj,tj,b,H,true)
        
        # compute probability density of going from merge launch state to original state
        log_prop_ba = update_parameter!(tm,ti0,xs,dummy,0,H,false)
        
        # compute acceptance probability
        log_prior_b = log_v[t+1] + lgamma(ni+b)+lgamma(nj+b)-2*lgamma(a) + log_prior(ti,H) + log_prior(tj,H)
        log_prior_a = log_v[t] + lgamma(ns+b)-lgamma(a) + log_prior(ti0,H)
        log_lik_ratio = 0.
        for k = 1:ns
            log_lik_ratio += log(likelihood(xs[k],(zs[k]? ti:tj))) - log(likelihood(xs[k],ti0))
        end
        p_accept = min(1, exp(log_prop_ba-log_prop_ab + log_prior_b-log_prior_a + log_lik_ratio))
        #println("split proposal: ",p_accept)
        
        # accept or reject
        accept = (rand()<p_accept)
        return accept,zs,ni,nj
        
    else  # propose a merge
        # make one final sweep and compute its probability density
        log_prop_ab = update_parameter!(tm,tm,xs,dummy,0,H,true)
        
        # compute probability density of going from split launch state to original state
        log_prop_ba,ni,nj = restricted_Gibbs!(zs,zs0,ni,nj,xs,ns,is,js,ti,ti0,tj,tj0,b,H,false)

        # compute acceptance probability
        log_prior_b = log_v[t-1] + lgamma(ns+b)-lgamma(a) + log_prior(tm,H)
        log_prior_a = log_v[t] + lgamma(ni+b)+lgamma(nj+b)-2*lgamma(a) + log_prior(ti0,H) + log_prior(tj0,H)
        @assert(ni==sum(zs0) && nj==sum(!zs0))
        log_lik_ratio = 0.
        for k = 1:ns
            log_lik_ratio += log(likelihood(xs[k],tm)) - log(likelihood(xs[k],(zs0[k]? ti0:tj0)))
        end
        p_accept = min(1, exp(log_prop_ba-log_prop_ab + log_prior_b-log_prior_a + log_lik_ratio))
        #println("merge proposal: ",p_accept)
        
        # accept or reject
        accept = (rand()<p_accept)
        return accept,zs,ns,ns
    end
end


function sampler(options,n_total,n_keep)
    x,n = options.x,options.n
    t_max = options.t_max
    a,b,log_v = options.a,options.b,options.log_v
    use_splitmerge,n_split,n_merge = options.use_splitmerge,options.n_split,options.n_merge
    use_hyperprior = options.use_hyperprior

    model_type = options.model_type
    alpha_random,p_alpha,alpha = options.alpha_random,options.p_alpha,options.alpha
    p_a = eval(parse(p_alpha))
    sigma_alpha = 0.1 # scale for MH proposals in alpha move

    @assert(n==length(x))
    keepers = zeros(Int,n_keep)
    keepers[:] = round(Int,linspace(round(Int,n_total/n_keep),n_total,n_keep))
    keep_index = 0

    if model_type=="MFM"
        A = a*exp(diff(log_v))
    else # DPM
        A = alpha*ones(n)
    end

    t = 1  # number of clusters
    z = ones(Int,n)  # z[i] = the cluster ID associated with data point i
    list = zeros(Int,t_max+3); list[1] = 1  # list[1:t] = the list of active cluster IDs
                                            # list is maintained in increasing order for 1:t, and is 0 after that.
                                            # Further, maximum(list[1:t]) is always <= t_max+3
    N = zeros(Int,t_max+3); N[1] = n  # N[c] = size of cluster c

    H = construct_hyperparameters(options)
    D = 2 # dimension of the parameter, theta = [mu,sigma]
    theta = zeros(D,t_max+3); theta[:,1] = prior_sample(H)  # theta[:,c] = parameters for cluster c
    theta_p = zeros(D)  # temporary variable to hold "proposed" theta
    theta_c = zeros(D)
    theta_aux = zeros(D,n)
    ti,tj,tm = zeros(D),zeros(D),zeros(D)

    indices = collect(1:n)  # used for sampling without replacement
    p = zeros(n+1)
    zs = ones(Int,n)  # temporary variable used for split-merge assignments
    S = zeros(Int,n)  # temporary variable used for split-merge indices
    
    log_Nb = log((1:n) + b)
    
    # Record-keeping variables
    t_r = zeros(Int8,n_total); @assert(t_max < 2^7)
    N_r = zeros(Int16,t_max+3,n_total); @assert(n < 2^15)
    z_r = zeros(Int8,n,n_keep); @assert(t_max < 2^7)
    theta_r = Theta[zeros(2) for j=1:t_max+3, iteration=1:n_keep]
    
    for iteration = 1:n_total
        #  -------------- Resample thetas and H --------------
        for j = 1:t
            c = list[j]
            theta_c[1] = theta[1,c]
            theta_c[2] = theta[2,c]
            update_parameter!(theta_c,theta_c,x,z,c,H,true)
            theta[1,c] = theta_c[1]
            theta[2,c] = theta_c[2]
        end
        if use_hyperprior
            update_hyperparameters!(H,theta,list,t)
        end
        if model_type=="DPM" && alpha_random
            # Metropolis-Hastings move for DP concentration parameter
            aprop = alpha*exp(randn()*sigma_alpha)
            top = t*log(aprop) - lgamma(aprop+n) + lgamma(aprop) + log(p_a(aprop)) + log(aprop)
            bot = t*log(alpha) - lgamma(alpha+n) + lgamma(alpha) + log(p_a(alpha)) + log(alpha)
            if rand() < min(1.0,exp(top-bot))
                alpha = aprop
            end
            #for i=1:t_max+1; log_v[i] = i*log(alpha) - lgamma(alpha+n) + lgamma(alpha); end
            log_v = (1:t_max+1)*log(alpha) - lgamma(alpha+n) + lgamma(alpha)
            A = alpha*ones(n)
        end
            
        # -------------- Resample z's --------------
        prior_sample!(theta_aux,n,H)
        for i = 1:n
            # remove point i from it's cluster
            c = z[i]    
            N[c] -= 1
            if N[c]>0
                # theta_p = theta_aux[:,i]
                theta_p[1] = theta_aux[1,i]
                theta_p[2] = theta_aux[2,i]
            else
                # theta_p = theta[:,c]
                theta_p[1] = theta[1,c]
                theta_p[2] = theta[2,c]
                # remove cluster {i}, keeping list in proper order
                ordered_remove!(c,list,t)
                t -= 1
            end
            
            # compute probabilities for resampling
            for j = 1:t
                p[j] = (N[list[j]]+b)*likelihood(x[i],theta,list[j])
            end
            p[t+1] = A[t]*likelihood(x[i],theta_p)
            
            # sample a new cluster for it
            j = randp(p,t+1)
            
            # add point i to it's new cluster
            if j<=t
                c = list[j]
            else
                @assert(t<t_max, "Sampled t exceeded t_max. Increase t_max and retry.")
                
                # give it the smallest available cluster ID number, and keep list ordered
                c = ordered_insert_next!(list,t)
                
                #theta[:,c] = theta_p
                theta[1,c] = theta_p[1]
                theta[2,c] = theta_p[2]
                t += 1
            end
            z[i] = c
            N[c] += 1
        end
        
        # -------------- Split/merge move --------------
        if use_splitmerge
            # randomly choose a pair of indices
            shuffle!(indices)
            i,j = indices[1],indices[2]
            c_i,c_j = z[i],z[j]
            # consider the subset of points in the clusters of i and j
            S = find((z.==c_i)|(z.==c_j))
            is,js = findfirst(S,i),findfirst(S,j)
            accept,zs,ni,nj = split_merge(x[S],(z[S].==c_i),is,js,ti,tj,tm,theta[:,c_i],theta[:,c_j],t,H,a,b,log_v,n_split,n_merge)
            if accept
                if c_i==c_j
                    # split off a new cluster, c_j
                    @assert(t<t_max, "Sampled t exceeded t_max. Increase t_max and retry.")
                    
                    # give the new cluster the smallest available ID number, and keep list ordered
                    c_j = ordered_insert_next!(list,t)
                    
                    # z[S[zs]] = c_i is already true
                    z[S[!zs]] = c_j
                    N[c_i] = ni
                    N[c_j] = nj
                    for l = 1:D
                        theta[l,c_i] = ti[l]
                        theta[l,c_j] = tj[l]
                    end
                    
                    t += 1
                else
                    # merge cluster c_j into c_i
                    z[S] = c_i
                    N[c_i] += N[c_j]
                    N[c_j] = 0  # (strictly speaking, this is not necessary)
                    theta[1,c_i] = tm[1]
                    theta[2,c_i] = tm[2]
                    # theta[:,c_j] = zeros(D)  # (not necessary)

                    # remove c_j from list, keeping list in proper order
                    ordered_remove!(c_j,list,t)
                    t -= 1
                end
            end
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
                theta_r[c,keep_index][1] = theta[1,c]
                theta_r[c,keep_index][2] = theta[2,c]
            end
        end
    end
    
    return t_r,N_r,z_r,theta_r,keepers
end

















