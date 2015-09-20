tmlenet
==========


[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tmlenet)](http://cran.r-project.org/package=tmlenet)
<!-- [![](http://cranlogs.r-pkg.org/badges/tmlenet)](http://cran.rstudio.com/web/packages/tmlenet/index.html) -->
[![Travis-CI Build Status](https://travis-ci.org/osofr/tmlenet.svg?branch=master)](https://travis-ci.org/osofr/tmlenet)
[![Coverage Status](https://coveralls.io/repos/osofr/tmlenet/badge.png?branch=master&service=github)](https://coveralls.io/r/osofr/tmlenet?branch=master)

The `tmlenet` R package implements the Targeted Maximum Likelihood Estimation (TMLE) of causal effects under single time point stochastic interventions in network data. The package also implements the Horvitz-Thompson estimator for networks and the parametric g-computation formula-based estimator. The inference for the TMLE is based on the efficient influence curves for dependent data. See the paper below for more information on the estimation methodology employed by `tmlenet`:

> M. J. van der Laan, “Causal inference for a population of causally connected units,” J. Causal Inference J. Causal Infer., vol. 2, no. 1, pp. 13–74, 2014.

### Installation

To install the development version of `tmlenet` (requires the `devtools` package):

```R
devtools::install_github('osofr/tmlenet', build_vignettes = FALSE)
```

### Documentation

Once the package is installed, please refer to the help file `?'tmlenet-package'` and `tmlenet` function documentation for details and examples:

```R
?'tmlenet-package'
?tmlenet
```

### Details

The input data are assumed to consist of unit-specific observations `(F,W,A,Y)`, where for each unit `i`: `F[i]` is a vector of the user-specified "__friends__" of unit `i` (also referred to as `i`'s __network__); `W[i]` are `i`'s baseline covariates; `A[i]` is `i`'s exposure (can be binary or continuous); and `Y[i]` is `i`'s outcome.

Each exposure `A[i]` can depend on baseline `W[i]` and the baseline covariate values of `i`'s friends, namely, all `W[j]` such that `j` is in set `F[i]`. Additionally, the outcome `Y[i]` can depend on `(W[i],A[i])` as well as the covariate values and exposures of `i`'s friends, i.e., all `(W[j], A[j])` such that `j` is in set `F[i]`.

The main function of the package is `tmlenet`. The input data format is a `data.frame` and it is assumed to consist of rows of unit-specific baseline covariates, exposures and outcomes. The network of friends can be specified either as a separate column in the input data (each `F[i]` is a string of space delimited IDs or row numbers of `i` friends) or as a separate matrix input argument to `tmlenet` function. Specifying the network via a separate input matrix generally leads to a significant improvement in run times. See `tmlenet` function help file for additional details on how to specifying these and the rest of the input arguments (such as, the stochastic interventions of interest, etc).


### Example
...

### Citation
To cite `tmlenet` in publications, please use:
> Sofrygin O, van der Laan MJ (2015). *tmlenet: Targeted Maximum Likelihood Estimation for Networks.* R package version 0.1.

### Funding
The development of this package was funded through an NIH grant (R01 AI074345-07).

### Copyright
This software is distributed under the GPL-2 license.
