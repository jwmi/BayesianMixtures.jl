# Gamma random number generator for Julia
# Original code by John D. Cook
# http://www.johndcook.com/julia_rng.html

module Gamma

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

end




