This folder contains examples using the BayesianMixtures package.
In particular, code is provided for the results in Miller and Harrison (2015).
Peter Green's Nmix program is also included.

----------------------------------------------------------------------
## Contents

#### Basic usage example
- basicExample.jl - A basic example of using the BayesianMixtures package.

#### Generating the results in Miller and Harrison (2015)
- galaxyCompareJNtoRJ.jl - (Section 7.1.1 of M&H) Compare Jain-Neal to reversible jump on galaxy data set.
- simCompareJNtoRJ.jl - (Section 7.1.2 of M&H) Compare Jain-Neal to reversible jump as dimensionality increases.
- geneExpression.jl - (Section 7.2 of M&H) Run Jain-Neal sampler on gene expression data.
- simCompareDPM.jl - (Section 7.3 of M&H) Compare MFM with DPM.
- datasets/ - Input data for these examples.

#### Green's Nmix program for reversible jump
- Nmix-source/ - Source code of Peter Green's Nmix program for RJMCMC.
- Nmix - Pre-compiled binary for Mac OS X. (NOTE: You will still need the gfortran libraries. See Nmix-source/README.md.)
- Nmix.exe - Pre-compiled executable for Windows. If you're lucky, this will work out-of-the-box.
- galx.dat - Galaxy data set in Nmix input format.
- galx/ - Folder to contain results from running Nmix on galx.dat.

----------------------------------------------------------------------


