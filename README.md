tmlenet
==========


<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tmlenet)](http://cran.r-project.org/package=tmlenet) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/tmlenet)](http://cran.rstudio.com/web/packages/tmlenet/index.html) -->
[![Travis-CI Build Status](https://travis-ci.org/osofr/tmlenet.svg?branch=master)](https://travis-ci.org/osofr/tmlenet)
[![Coverage Status](https://coveralls.io/repos/osofr/tmlenet/badge.png?branch=master&service=github)](https://coveralls.io/r/osofr/tmlenet?branch=master)

The `tmlenet` R package performs estimation of average causal effects for single time point interventions in network-dependent (non-IID) data in the presence of interference and/or spillover. Currently implemented estimation algorithms are the targeted maximum likelihood estimation (TMLE), Horvitz-Thompson or the inverse-probability-of-treatment (IPTW) estimator and the parametric G-computation estimator. The user-specified interventions can be either static, dynamic or stochastic. Asymptotically correct influence-curve-based confidence intervals are also constructed for the TMLE and IPTW. See the paper below for more information on the estimation methodology employed by the `tmlenet` R package:

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

### The input data and the network summary measures

The input data are assumed to consist of rows of unit-specific observations, with each row `i` represented by variables (`F.i`,`W.i`,`A.i`,`Y.i`), where `F.i` is a vector of "__friend IDs__" of unit `i` (also referred to as `i`'s "__network__"), `W.i` is a vector of `i`'s baseline covariates, `A.i` is `i`'s exposure (either binary, categorical or continuous) and `Y.i` is `i`'s binary outcome. 

Each exposure `A.i` depends on (possibly multivariate) baseline summary measure(s) `sW.i`, where `sW.i` can be any user-specified function of `i`'s baseline covariates `W.i` and the baseline covariates of `i`'s friends in `F.i` (all `W.j` such that `j` is in `F.i`). Similarly, each outcome `Y.i` depends on `sW.i` and (possibly multivariate) summary measure(s) `sA.i`, where `sA.i` can be any user-specified function of `i`'s baseline covariates and exposure (`W.i`,`A.i`) and the baseline covariates and exposures of `i`'s friends (all `W.j`,`A.j` such that `j` is in `i`'s friend set `F.i`). 

The summary measures (`sW.i`,`sA.i`) are defined simultaneously for all `i` with functions `def.sW` and `def.sA`. It is assumed that (`sW.i`,`sA.i`) have the same dimensionality across `i`. The function `eval.summaries` can be used for evaluating these summary measures.

All estimation is performed by calling the `tmlenet` function. The vector of friends `F.i` can be specified either as a single column in the input data (where each `F.i` is a string of friend IDs or friend row numbers delimited by character `sep`) or as a separate input matrix of network IDs (where each row is a vector of friend IDs or friend row numbers). Specifying the network as a matrix generally results in significant improvements to run time. See `tmlenet` function help file for additional details on how to specify these and the rest of the input arguments.

### Example
...

### Citation
To cite `tmlenet` in publications, please use:
> Sofrygin O, van der Laan MJ (2015). *tmlenet: Targeted Maximum Likelihood Estimation for Networks.* R package version 0.1.

### Funding
The development of this package was funded through an NIH grant (R01 AI074345-07).

### Copyright
This software is distributed under the GPL-2 license.
