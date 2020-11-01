# Random number generators for Julia
# This is a slightly modified and reduced version of a Julia module by John D. Cook.
# http://www.johndcook.com/julia_rng.html

module RandomNumbers

inverse_gamma(shape, rate) = 1/gamma(shape, rate)

function gamma(shape, rate)
    # Return a random sample from a gamma distribution
    # with pdf \propto x^{shape-1}exp(-rate*x)

    if shape <= 0.0
        error("Shape parameter must be positive, but shape = $shape")
    end
    if rate <= 0.0
        error("Rate parameter must be positive, but rate = $rate")
    end
    
    ## Implementation based on "A Simple Method for Generating Gamma Variables"
    ## by George Marsaglia and Wai Wan Tsang.  
    ## ACM Transactions on Mathematical Software
    ## Vol 26, No 3, September 2000, pages 363-372.

    if shape >= 1.0
        d = shape - 1.0/3.0
        c = 1.0/sqrt(9.0*d)
        while true
            x = randn()
            v = 1.0 + c*x
            while v <= 0.0
                x = randn()
                v = 1.0 + c*x
            end
            v = v*v*v
            u = rand()
            xsq = x*x
            if u < 1.0 -.0331*xsq*xsq || log(u) < 0.5*xsq + d*(1.0 - v + log(v))
                return d*v/rate
            end
        end
    else
        g = gamma(shape+1.0, 1.0)
        w = rand()
        return g*(w^(1.0/shape))/rate
    end
end

## return a random sample from a chi square distribution
## with the specified degrees of freedom
function chi_square(dof)
    return gamma(0.5*dof, 0.5)
    # This was corrected from gamma(0.5, 2.0*dof) in the original code
end

## return a random sample from a beta distribution
function beta(a, b)
    if (a <= 0) || (b <= 0)
        error("Beta parameters must be positive")
    end
    
    ## There are more efficient methods for generating beta samples.
    ## However such methods are a little more efficient and much more complicated.
    ## For an explanation of why the following method works, see
    ## http://www.johndcook.com/distribution_chart.html#gamma_beta

    u = gamma(a, 1.0)
    v = gamma(b, 1.0)
    return u / (u + v)
end

end




