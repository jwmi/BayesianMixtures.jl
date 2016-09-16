# BayesianMixtures

[![Build Status](https://travis-ci.org/jwmi/BayesianMixtures.jl.svg?branch=master)](https://travis-ci.org/jwmi/BayesianMixtures.jl)

## About

BayesianMixtures is a Julia package for nonparametric Bayesian mixture models. The following model types are currently implemented:
- mixtures with a prior on the number of components, and
- Dirichlet process mixtures.

The following component distributions are currently implemented:
- univariate normal,
- multivariate normal with diagonal covariance matrix, and
- multivariate normal (with unconstrained covariance matrix).

For all models, inference is performed using the Jain-Neal split-merge sampler, including both conjugate and non-conjugate cases (Jain and Neal, 2004, 2007).  For mixtures with a prior on the number of components, this is done using the results of Miller and Harrison (2016).

Please cite the following publication if you use this package for your research:
> J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *arXiv preprint*, http://arxiv.org/abs/1502.06241, 2016.


## Installation

- Install [Julia](http://julialang.org/downloads/).

- Start Julia and run the following command:
```
Pkg.clone("https://github.com/jwmi/BayesianMixtures.jl.git")
```

### OPTIONAL

For full functionality, [JLD](https://github.com/JuliaIO/JLD.jl) and [PyPlot](https://github.com/JuliaPy/PyPlot.jl) are also recommended. JLD enables you to save/load results from file. PyPlot enables you to plot results. To install JLD, run the following command in Julia:
```
Pkg.add("JLD")
```
Installing PyPlot can be complicated, but hopefully the following commands will work:
```
ENV["PYTHON"]=""
Pkg.add("PyPlot")
```
If not, please see the [PyPlot installation instructions](https://github.com/JuliaPy/PyPlot.jl).


## Basic usage example

```julia
using BayesianMixtures

# Simulate some data
x = randn(500)

# Specify model, data, and MCMC options
n_total = 1000  # total number of MCMC sweeps to run
options = BayesianMixtures.options("Normal","MFM",x,n_total)  # MFM model with univariate Normal components

# Run MCMC sampler
result = BayesianMixtures.run_sampler(options)

# Get the posterior on k (number of components) 
k_posterior = BayesianMixtures.k_posterior(result)
```


## Additional features

In addition to the Jain-Neal sampler, reversible jump MCMC is also implemented for certain classes of mixtures with a prior on the number of components (univariate normal mixtures and multivariate normal mixtures with diagonal covariance). For univariate normal mixtures, we use the algorithm of Richardson and Green (1997). A copy of Peter Green's [Nmix](https://people.maths.bris.ac.uk/~mapjg/Nmix/) program is also included, for convenience.


## References

S. Jain and R. M. Neal. A split-merge Markov chain Monte Carlo procedure for the Dirichlet process mixture model. *Journal of Computational and Graphical Statistics*, 13(1), 2004.

S. Jain and R. M. Neal. Splitting and merging components of a nonconjugate Dirichlet process mixture model. *Bayesian Analysis*, 2(3):445-472, 2007.

J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *arXiv preprint*, http://arxiv.org/abs/1502.06241, 2016.

S. Richardson and P. J. Green. On Bayesian analysis of mixtures with an unknown number of components. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 59(4):731-792, 1997.



