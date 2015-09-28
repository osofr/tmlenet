tmlenet
==========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tmlenet)](http://cran.r-project.org/package=tmlenet)
[![](http://cranlogs.r-pkg.org/badges/tmlenet)](http://cran.rstudio.com/web/packages/tmlenet/index.html)
[![Travis-CI Build Status](https://travis-ci.org/osofr/tmlenet.svg?branch=master)](https://travis-ci.org/osofr/tmlenet)
[![Coverage Status](https://coveralls.io/repos/osofr/tmlenet/badge.png?branch=master&service=github)](https://coveralls.io/r/osofr/tmlenet?branch=master)

The `tmlenet` R package performs estimation of average causal effects for single time point interventions in network-dependent (non-IID) data in the presence of interference and/or spillover. Currently implemented estimation algorithms are the targeted maximum likelihood estimation (TMLE), Horvitz-Thompson or the inverse-probability-of-treatment (IPTW) estimator and the parametric G-computation estimator. The user-specified interventions can be either static, dynamic or stochastic. Asymptotically correct influence-curve-based confidence intervals are also constructed for the TMLE and IPTW. See the paper below for more information on the estimation methodology employed by the `tmlenet` R package:

> M. J. van der Laan, “Causal inference for a population of causally connected units,” J. Causal Inference J. Causal Infer., vol. 2, no. 1, pp. 13–74, 2014.

### Installation

To install the CRAN release version of `simcausal`: 

```R
install.packages('tmlenet')
```

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

We will use the sample dataset (`W`=(`W1`,`W2`,`W3`),`A`,`Y`) and the sample network matrix of friend IDs (`F`) that come along with the package:

```R
data(df_netKmax6)
head(df_netKmax6)
data(NetInd_mat_Kmax6)
head(NetInd_mat_Kmax6)
Kmax <- ncol(NetInd_mat_Kmax6) # Max number of friends in this network:
```

The estimation algorithm assumes that the outcomes in `Y.i` for units `i=1,...,N` are conditionally independent,
given the summary measures defined in `def_sW` and the summary measures defined in `def_sA`.

When no additional assumptions about the conditional independence of outcomes `Y.i` can be made 
(beyond the dependence on the network structure),
one can define the summary measures `sW` and `sA` non-parametrically, e.g., 
for each observation `i`: include in `sW` all baseline covariates of unit `i` and 
all baseline covariates of `i`'s friends; include in `sA` the exposure of unit `i` and
all exposures of `i`'s friends. 

The example below does just that, defining  `sW`:=(`netW1`,`netW2`,`netW3`) and `sA`:=`netA`, 
where `netVar` is a summary measure of dimension `Kmax+1` and includes `Var` values of each 
unit as well as `Var` values of all friends of each unit:

```R
def_sW <- def.sW(netW1 = W1[[0:Kmax]], netW2 = W2[[0:Kmax]], netW3 = W3[[0:Kmax]])
def_sA <- def.sA(netA = A[[0:Kmax]])
```

Note that the summary measure `nF` (number of friends for each unit) is always added automatically to 
`def.sW` function calls (only once), but not to `def.sA`.

A helper function that can pre-evaluate the above summary measures based on the input data:

```R
eval_res <- eval.summaries(sW = def_sW, sA = def_sA,  Kmax = 6, data = df_netKmax6,
                          NETIDmat = NetInd_mat_Kmax6)
```

Contents of the list returned by eval.summaries():

```R
head(eval_res$sW.matrix) # Matrix of sW summary measures:
head(eval_res$sA.matrix) # Matrix of sA summary measures:
head(eval_res$NETIDmat) # matrix of network IDs:
# Observed data summary measures (sW,sA) and network stored in one object:
# eval_res$DatNet.ObsP0
# class(eval_res$DatNet.ObsP0)
```

In the example below, we estimate mean population outcome under deterministic intervention that assigns all `A` to 0
(network specified via a matrix of friend IDs). Note that can also use previously evaluated
summary measures object `DatNet.ObsP0` as input to `tmlenet`, avoiding the need to specify the arguments
(`data`,`NETIDmat`,`Kmax`,`sW`,`sA`) for the second time.

```R
res1 <- tmlenet(data = df_netKmax6, NETIDmat = NetInd_mat_Kmax6, Kmax = Kmax, 
                sW = def_sW, sA = def_sA,
                Anode = "A", Ynode = "Y",
                f_gstar1 = 0L, optPars = list(n_MCsims = 1))
res1$EY_gstar1$estimates
res1$EY_gstar1$vars
res1$EY_gstar1$CIs
```

By default, the conditional expectation `E[Y=1|...]` (`Qform` argument) is estimated by including all
summary measures in `sW` and `sA` as predictors in the logistic regression for the outcome `Y`.
Similarly, by default, the observed exposure model `P(sA|sW)` (`hform.g0` argument) is estimated
as the conditional probability of observing the summary measures defined in `sA`, given the summary measures
defined in `sW`. Finally, the intervention exposure model `P(sA^*|sW)` (`hform.gstar` 
argument) is estimated by first replacing all observed exposures in `A` with those generated from
the intervention function specified in `f_gstar1` (new exposures denoted by `A^*`) and then building
the same summary measures defined in `sW` and `sA` using exposures `A^*` instead of `A`
(new summary measures denoted by `sA^*`). By default, the intervention exposure model `P(sA^*|sW)`
will be estimated as the conditional probability of observing the intervention-based summary measures in `sA^*` 
(`sA^*` built with `A^*` using the same summary mappings as in `sA`), given the summary measures defined in `sW`.

One can change this default behavior and use the arguments `Qform`, `hform.g0` and `hform.gstar`
to select a subset of the summary measures in `sW`,`sA` to be included in each of the three models described above.
For example, below we are assuming that the outcomes in `Y` only depend on the summary measures `netA`,`netW2` 
(regression `"Y~netA+netW2"`) and hence the observed exposure model is given by `P(netA|netW2)` (regression `"netA~netW2"`) and we also know that `f_gstar1` defines a static intervention `A^*=1` and hence `sA^*` is degenerate and doesn't depend on any baseline covariates and will be estimated here with a simplified regression model (regression `"netA ~ nF"`):

```R
res2 <- tmlenet(DatNet.ObsP0 = eval_res$DatNet.ObsP0,
                    Anode = "A", Ynode = "Y", 
                    Qform = "Y ~ netA + netW2",
                    hform.g0 = "netA ~ netW2",
                    hform.gstar = "netA ~ nF",
                    f_gstar1 = 0L, optPars = list(n_MCsims = 1))
res2$EY_gstar1$estimates
res2$EY_gstar1$vars
res2$EY_gstar1$CIs
```

One might be also willing to make dimension reducing assumptions about the dependence of each `Y.i` on its 
network. For example, here we assume that each `Y.i` depends on its network's baseline covariates only
through a sum of its friends' values of `W3` and `Y.i` depends on its network's exposures only through a sum 
of `i`'s friends' interactions `(1-A)*(W2)` (while we assume `Y.i` still depends on `i`'s baseline covariates and
`i`'s exposure):

```R
def_sW <- def.sW(W = c(W1,W2,W3)) +
          def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

def_sA <- def.sA(A) +
          def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE)

eval_res <- eval.summaries(sW = def_sW, sA = def_sA, Kmax = 6, data = df_netKmax6,
                            NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

res3 <- tmlenet(DatNet.ObsP0 = eval_res$DatNet.ObsP0,
                Anode = "A", Ynode = "Y",
                Qform = "Y ~ A + sum.netAW2 + W + sum.netW3 + nF",
                hform.g0 = "A + sum.netAW2 ~ sum.netW3",
                hform.gstar = "A + sum.netAW2 ~ sum.netW3",
                f_gstar1 = 0, optPars = list(n_MCsims = 1))
res3$EY_gstar1$estimates
```

Note that the above model specified by `Qform` includes all summary measures in `sW`,`sA`, and hence is equivalent to the default regression model that would have been used if `Qform` was omitted.

One can specify any intervention of interest, for example below we estimate the counterfactual mean outcome under intervention that randomly assigns 20% of the population to exposure `A=1`. Note that we are also increasing the number of Monte-Carlo simulations
from 1 to 100.

```R
f.A_.2 <- function(data, ...) rbinom(n = nrow(data), size = 1, prob = 0.2)
res4 <- tmlenet(data = df_netKmax6, NETIDmat = NetInd_mat_Kmax6, Kmax = Kmax,
                sW = def_sW, sA = def_sA, 
                Anode = "A", Ynode = "Y", 
                f_gstar1 = f.A_.2, optPars = list(n_MCsims = 100))
res4$EY_gstar1$estimates
```

To estimate the average treatment effect (ATE) for two interventions (static or stochastic), specify the second intervention function using the argument `optPars(f_gstar2 = ...)`. In the example below, the intervention `f_gstar1`
statically sets everyone's exposure to `A=1` and the intervention `f_gstar2` statically sets everyone's exposure to `A=0`:

```R
res5 <- tmlenet(data = df_netKmax6, NETIDmat = NetInd_mat_Kmax6, Kmax = Kmax,
                sW = def_sW, sA = def_sA, Anode = "A", Ynode = "Y",
                f_gstar1 = 1, optPars = list(f_gstar2 = 0, n_MCsims = 1))
res5$ATE$estimates
```

### Citation
To cite `tmlenet` in publications, please use:
> Sofrygin O, van der Laan MJ (2015). *tmlenet: Targeted Maximum Likelihood Estimation for Networks.* R package version 0.1.

### Funding
The development of this package was funded through an NIH grant (R01 AI074345-07).

### Copyright
This software is distributed under the GPL-2 license.
