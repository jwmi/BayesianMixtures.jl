# BayesianMixtures

<!--
[![Build Status](https://travis-ci.org/jwmi/BayesianMixtures.jl.svg?branch=master)](https://travis-ci.org/jwmi/BayesianMixtures.jl)
-->

## About

BayesianMixtures is a Julia package for nonparametric Bayesian mixture models and Bayesian clustering. The following model types are currently implemented:
- mixture with a prior on the number of components, a.k.a., mixture of finite mixtures (MFM), and
- Dirichlet process mixture.

The following component distributions are currently implemented:
- univariate normal,
- multivariate normal with diagonal covariance matrix, and
- multivariate normal (with unconstrained covariance matrix).

For all models, inference is performed using the Jain-Neal split-merge samplers, including both conjugate and non-conjugate cases (Jain and Neal, 2004, 2007).  For MFMs, this is done using the results of Miller and Harrison (2017).

Please cite the following publication if you use this package in your research:
> J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *Journal of the American Statistical Association*, 2017. (In press) http://arxiv.org/abs/1502.06241.


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

### Optional: Saving results

Functions for saving/loading results from file are provided.  To use them, you need to install a couple things as follows:
- Mac users (Windows users skip this step): Install the Command Line Tools (CLT) by opening a terminal window and entering: `xcode-select --install`.  In the window that pops up, click `Install` (**not** `Get Xcode`).
- Then, install [JLD](https://github.com/JuliaIO/JLD.jl) by opening Julia and entering: `Pkg.add("JLD")`.

### Optional: Plotting results

Functions for plotting are also provided.  To use them, install [PyPlot](https://github.com/JuliaPy/PyPlot.jl) by entering the following commands:
```
ENV["PYTHON"]=""
if (Pkg.installed("PyCall")!=nothing); Pkg.build("PyCall"); end
Pkg.add("PyPlot")
using PyPlot
```
Hopefully that will work, but if you get errors, you can either (a) not use the plotting functions in BayesianMixtures, or (b) try your luck with the PyPlot installation instructions [here](https://github.com/JuliaPy/PyPlot.jl).  



## Updating to the latest version

If you want to update your installation of BayesianMixtures to the latest version of the package, run the following commands in Julia:
```
Pkg.rm("BayesianMixtures")
Pkg.clone("https://github.com/jwmi/BayesianMixtures.jl.git")
```


## Questions or bugs

If you have a question or find a bug, feel free to contact me ([Jeff Miller](http://jwmi.github.io/)). Also feel free to submit a pull request if you find and fix a bug.


## Licensing / Citation

This package is released under an MIT license (with the exception of Peter Green's Nmix code and John D. Cook's random number generators). See [LICENSE.md](LICENSE.md). 

Please cite the following publication if you use this package in your research:
J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *Journal of the American Statistical Association*, 2017. (In press) http://arxiv.org/abs/1502.06241.


## References

S. Jain and R. M. Neal. A split-merge Markov chain Monte Carlo procedure for the Dirichlet process mixture model. *Journal of Computational and Graphical Statistics*, 13(1), 2004.

S. Jain and R. M. Neal. Splitting and merging components of a nonconjugate Dirichlet process mixture model. *Bayesian Analysis*, 2(3):445-472, 2007.

J. W. Miller and M. T. Harrison. Mixture models with a prior on the number of components. *Journal of the American Statistical Association*, 2017. (In press) http://arxiv.org/abs/1502.06241.

S. Richardson and P. J. Green. On Bayesian analysis of mixtures with an unknown number of components. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 59(4):731-792, 1997.

