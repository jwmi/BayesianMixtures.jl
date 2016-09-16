This folder contains the source code for the BayesianMixtures package.

----------------------------------------------------------------------
## Contents

#### Main module file
- BayesianMixtures.jl - Specify options, call sampler, analyze results, plot results.

#### Samplers
- coreConjugate.jl - collapsed Jain-Neal sampler for conjugate priors.
- coreNonconjugate.jl - (uncollapsed) Jain-Neal sampler for non-conjugate priors.
- coreNormal.jl - (uncollapsed) Jain-Neal sampler optimized for univariate normal model.

#### Component distributions
- Normal.jl - univariate normal with parameter settings from Richardson and Green (1997), using coreNormal.jl sampler.
- MVN.jl - multivariate normal (with unconstrained covariance matrices), using coreNonconjugate.jl.
- MVNaaC.jl - multivariate normal with diagonal covariance matrices, using coreConjugate.jl.
- MVNaaN.jl - multivariate normal with diagonal covariance matrices, using coreNonconjugate.jl.
- MVNaaRJ.jl - multivariate normal with diagonal covariance matrices, using reversible jump sampler.
- NormalNonoptimized.jl - univariate normal with parameter settings from Richardson and Green (1997), using generic coreNonconjugate.jl.

#### Partition distributions
- MFM.jl - mixture of finite mixtures (MFM) model.

#### Helper functions
- Gamma.jl - sampling gamma random variables (faster than calling Distributions.jl).
- Lower.jl - matrix operations using lower-triangular matrices (faster than built-in Julia matrix operations).



----------------------------------------------------------------------


