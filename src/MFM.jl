# Functions for computing MFM partition distribution coefficients V_n(t) and also p(k|t).
module MFM

logsumexp(a,b) = (m = max(a,b); m == -Inf ? -Inf : log(exp(a-m) + exp(b-m)) + m)

# Compute log_v[t] = log(V_n(t)) under the given MFM parameters, for t=1:upto.
function coefficients(log_pk,gamma,n,upto)
    tolerance = 1e-12
    log_v = zeros(upto)
    for t = 1:upto
        if t>n; log_v[t] = -Inf; continue; end
        a,c,k,p = 0.0, -Inf, 1, 0.0
        while abs(a-c) > tolerance || p < 1.0-tolerance  # Note: The first condition is false when a = c = -Inf
            if k >= t
                a = c
                b = lgamma(k+1)-lgamma(k-t+1) - lgamma(k*gamma+n)+lgamma(k*gamma) + log_pk(k)
                c = logsumexp(a,b)
            end
            p += exp(log_pk(k))
            k = k+1
        end
        log_v[t] = c
    end
    return log_v
end

# Compute p[k,t]=p(k|t) under the given MFM parameters, for k=1:upto_k, t=1:upto_t.
function p_kt(log_pk,gamma,n,upto_k,upto_t)
    @assert(upto_t<=n,"p(k|t) is undefined for t>n.")
    log_v = coefficients(log_pk,gamma,n,upto_t)
    p = zeros(upto_k,upto_t)
    for k = 1:upto_k, t = 1:min(k,upto_t)
        b = lgamma(k+1)-lgamma(k-t+1) - lgamma(k*gamma+n)+lgamma(k*gamma) + log_pk(k)
        p[k,t] = exp(b - log_v[t])
    end
    return p
end


end # module MFM

