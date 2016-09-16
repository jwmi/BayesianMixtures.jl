# BayesianMixtures

<!--
[![Build Status](https://travis-ci.org/jwmi/BayesianMixtures.jl.svg?branch=master)](https://travis-ci.org/jwmi/BayesianMixtures.jl)
-->

## About

BayesianMixtures is a Julia package for nonparametric Bayesian mixture models and Bayesian clustering. The following model types are currently implemented:
- mixtures with a prior on the number of components, a.k.a., mixture of finite mixtures (MFM), and
- Dirichlet process mixtures.

The following component distributions are currently implemented:
- univariate normal,
- multivariate normal with diagonal covariance matrix, and
- multivariate normal (with unconstrained covariance matrix).

For all models, inference is performed using the Jain-Neal split-merge sampler, including both conjugate and non-conjugate cases (Jain and Neal, 2004, 2007).  For mixtures with a prior on the number of components, this is done using the results of Miller and Harrison (2015).

Please cite the following publication if you use this package for your research:
> J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *arXiv preprint*, http://arxiv.org/abs/1502.06241, 2015.


## Installation

- Install [Julia](http://julialang.org/downloads/).

- Start Julia and run the following command:
```
Pkg.clone("https://github.com/jwmi/BayesianMixtures.jl.git")
```

## Basic usage example

```julia
using BayesianMixtures

# Simulate some data
x = randn(500)

# Specify model, data, and MCMC options
n_total = 1000  # total number of MCMC sweeps to run
options = BayesianMixtures.options("Normal","MFM",x,n_total)  # MFM with univariate Normal components

# Run MCMC sampler
result = BayesianMixtures.run_sampler(options)

# Get the posterior on k (number of components) 
k_posterior = BayesianMixtures.k_posterior(result)
```

For more in-depth examples, see the [examples](examples/) folder.


## Additional features

### Optional: Reversible jump MCMC

In addition to the Jain-Neal sampler, reversible jump MCMC is also implemented for certain classes of MFMs (specifically, univariate normal mixtures and multivariate normal mixtures with diagonal covariance). For univariate normal mixtures, we use the algorithm of Richardson and Green (1997). A copy of Peter Green's [Nmix](https://people.maths.bris.ac.uk/~mapjg/Nmix/) program is also included, for convenience.

### Optional: Plotting results

Functions for plotting are also provided.  To use them, you need to install [PyPlot](https://github.com/JuliaPy/PyPlot.jl).  Installing PyPlot can be complicated, but hopefully the following commands will work:
```
ENV["PYTHON"]=""
Pkg.add("PyPlot")
```
If that doesn't work, you can either (a) not use the plotting functions in BayesianMixtures, or (b) try your luck with the [PyPlot installation instructions](https://github.com/JuliaPy/PyPlot.jl).


## Updating to the latest version

If you want to update your installation of BayesianMixtures to the latest version of the package, run the following commands in Julia:
```
Pkg.rm("BayesianMixtures")
Pkg.clone("https://github.com/jwmi/BayesianMixtures.jl.git")
```


## Contact / bug fixes

If you have a question or find a bug, feel free to contact me ([Jeff Miller](http://jwmi.github.io/)). Also feel free to submit a pull request if you find and fix a bug.


## References

S. Jain and R. M. Neal. A split-merge Markov chain Monte Carlo procedure for the Dirichlet process mixture model. *Journal of Computational and Graphical Statistics*, 13(1), 2004.

S. Jain and R. M. Neal. Splitting and merging components of a nonconjugate Dirichlet process mixture model. *Bayesian Analysis*, 2(3):445-472, 2007.

J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *arXiv preprint*, http://arxiv.org/abs/1502.06241, 2015.

S. Richardson and P. J. Green. On Bayesian analysis of mixtures with an unknown number of components. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 59(4):731-792, 1997.



