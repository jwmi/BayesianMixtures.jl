# Lower triangular matrices
module Lower

include("Gamma.jl")
using .Gamma

function Cholesky!(L,A,d)
    # Compute L lower triangular such that L*L' = A, given a d-by-d positive definite matrix A.
    # Only L is modified by this procedure.
    for i = 1:d
        s = 0.
        for j = 1:i-1
            r = 0.
            for k = 1:j-1; r += L[i,k]*L[j,k]; end
            L[i,j] = (A[i,j]-r)/L[j,j]
            s += L[i,j]*L[i,j]
        end
        L[i,i] = sqrt(A[i,i]-s)
        for j = i+1:d; L[i,j] = 0.; end
    end
end

function reverse_Cholesky!(L,A,d)
    # Compute L lower triangular such that L'*L = A, given a d-by-d positive definite matrix A.
    # Only L is modified by this procedure.
    for j = d:-1:1
        s = 0.
        for i = d:-1:j+1
            r = 0.
            for k = i+1:d; r += L[k,i]*L[k,j]; end
            L[i,j] = (A[i,j]-r)/L[i,i]
            s += L[i,j]*L[i,j]
        end
        L[j,j] = sqrt(A[j,j]-s)
        for i = 1:j-1; L[i,j] = 0.; end
    end
end

function Cholesky_inverse!(L,M,A,d)
    # Compute L and M lower triangular such that L*L' = inv(A) and M'*M = A,
    # given a d-by-d positive definite matrix A. Only L and M are modified by this procedure.
    reverse_Cholesky!(M,A,d) # M'*M = A
    inverse!(L,M,d) # L = inv(M)
end

function scale!(L,c,d)
    # Multiply L by c, for a scalar c and a d-by-d matrix L.
    for i = 1:d*d; L[i] *= c; end
end

function multiply!(L,M,N,d)
    # Compute L = M*N, for d-by-d lower triangular matrices M and N.
    # NOTE: This is carefully implemented so that L and N can occupy the same memory.
    for i = d:-1:1, j = i:-1:1
        r = 0.
        for k = j:i
            r += M[i,k]*N[k,j]
        end
        L[i,j] = r
    end
end

function multiplyMNt!(L,M,N,d)
    # Compute L = M*N', for d-by-d lower triangular matrices M and N.
    for i = 1:d, j = 1:d
        r = 0.
        for k = 1:min(i,j)
            r += M[i,k]*N[j,k]
        end
        L[i,j] = r
    end
end

function multiplyMtN!(L,M,N,d)
    # Compute L = M'*N, for d-by-d lower triangular matrices M and N.
    # NOTE: This is carefully implemented so that L and N can occupy the same memory.
    for i = 1:d
        for j = 1:d
            r = 0.
            for k = max(i,j):d
                r += M[k,i]*N[k,j]
            end
            L[i,j] = r
        end
    end
end

function logdetsq(L,d)
    # Return log(det(L*L'))=log(det(L'*L))=log(det(L)^2) for a d-by-d triangular matrix L.
    r = 0.
    for i = 1:d; r += log(L[i,i]); end
    return 2.0*r
end

function quadratic(x,y,L,d)
    # Return (x-y)'*L*L'*(x-y) for vectors x,y and a d-by-d lower triangular matrix L.
    r = 0.
    for j = 1:d
        s = 0.
        for i = j:d; s += (x[i]-y[i])*L[i,j]; end
        r += s*s
    end
    return r
end

function sample!(L,M,nu,d)
    # Fill L with a randomly generated lower triangular matrix such that L*L' ~ Wishart(M*M',nu),
    # that is, L*L' is Wishart with scale matrix V=M*M' and nu degrees of freedom.
    # L and M are d-by-d lower triangular matrices, and nu is a real number greater than d-1.
    # 
    # This uses the Bartlett decomposition.
    # The same approach is used by Matlab's wishrnd function, so it should be valid 
    # for any nu > d-1, including noninteger nu.
    @assert(nu>d-1)
    for i = 1:d
        for j = 1:i-1; L[i,j] = randn(); end
        L[i,i] = sqrt(Gamma.chi_square(nu-i+1))
        for j = i+1:d; L[i,j] = 0.; end
    end
    multiply!(L,M,L,d)
    return L
end

function inverse!(L,M,d)
    # Compute L = inv(M), for an invertible d-by-d lower triangular matrix M.
    # NOTE: This is faster than inv(M) for d up to around 100, but slower after that (as of Dec 18, 2014).
    for i = 1:d
        for j = i+1:d; L[i,j] = 0.; end
        @assert(M[i,i] != 0.)
        L[i,i] = 1.0/M[i,i]
        for j = i-1:-1:1
            r = 0.
            for k = j+1:i
                r += L[i,k]*M[k,j]
            end
            L[i,j] = -r*L[j,j]  # Multiplication is faster than division
        end
    end
end

function solve_Lx_eq_y!(L,x,y,d)
    # Solve L*x = y for x, for a d-by-d lower triangular matrix L. Only x is modified by this procedure.
    # This uses the forward-substitution algorithm, avoiding the need to invert L.
    # The algorithm is valid if x and y occupy the same memory.
    for i = 1:d
        s = 0.
        for j = 1:i-1
            s += L[i,j]*x[j]
        end
        x[i] = (y[i]-s)/L[i,i]
    end
end

function solve_Ltx_eq_y!(L,x,y,d)
    # Solve L'*x = y for x, for a d-by-d lower triangular matrix L. Only x is modified by this procedure.
    # This uses the back-substitution algorithm, avoiding the need to invert L'.
    # The algorithm is valid if x and y occupy the same memory.
    for i = d:-1:1
        s = 0.
        for j = i+1:d
            s += L[j,i]*x[j]
        end
        x[i] = (y[i]-s)/L[i,i]
    end
end

function sample_Normal!(x,m,L,d)
    # Sample x ~ Normal(m,inv(L*L')), where L is an invertible d-by-d lower triangle matrix.
    # This is done by solving L'*(x-m) = z for x, where z ~ N(0,I). Only x is modified by this procedure.
    for i = d:-1:1
        s = 0.
        for j = i+1:d
            s += L[j,i]*x[j]
        end
        x[i] = (randn()-s)/L[i,i]
    end
    for i = 1:d; x[i] += m[i]; end
    return x
end











end # module Lower


