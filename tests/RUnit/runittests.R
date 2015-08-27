### --- Test setup ---
`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))

 
if(FALSE) {
  library("RUnit")
  library("roxygen2")
  library("devtools")
  setwd(".."); setwd(".."); getwd()
  document()
  load_all("./") # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # tmlenet:::debug_set() # SET TO DEBUG MODE

  setwd("..");
  install("tmlenet", build_vignettes = FALSE) # INSTALL W/ devtools:

  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf tmlenet")  # just create the pdf manual from help files

  # CHECK AND BUILD PACKAGE:
  getwd()
  # setwd("./tmlenet"); setwd(".."); getwd()
  devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs check with devtools
  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  devtools::build(args = "--compact-vignettes") # build package tarball compacting vignettes
  # devtools::build(args = "--no-build-vignettes") # build package tarball compacting vignettes
  # devtools::build() # build package tarball
  setwd("..")
  system("R CMD check --as-cran tmlenet_0.2.0.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes tmlenet") # check without building the pdf manual and not building vignettes
      ## system("R CMD build tmlenet --no-build-vignettes")
      ## system("R CMD build tmlenet")  
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./tmlenet_0.2.0.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(tmlenet)
  # tmlenet:::debug_set() # SET TO DEBUG MODE
  # tmlenet:::debug_off() # SET DEBUG MODE OFF

  # To install a specific branch:
  # devtools::install_github('osofr/simcausal', ref = "simnet", build_vignettes = FALSE)
  # devtools::install_github('osofr/tmlenet', ref = "master", build_vignettes = FALSE)

  # TEST COVERATE:
  # if your working directory is in the packages base directory
  # package_coverage()
  # or a package in another directory
  # cov <- package_coverage("tmlenet")
  # view results as a data.frame
  # as.data.frame(cov)
  # zero_coverage() can be used to filter only uncovered lines.
  # zero_coverage(cov)
}

sample_checks <- function() {   # doesn't run, this is just to show what test functions can be used
  print("Starting tests...")
    checkTrue(1 < 2, "check1")     ## passes fine
     ## checkTrue(1 > 2, "check2")  ## appears as failure in the test protocol
     v <- 1:3
     w <- 1:3
     checkEquals(v, w)               ## passes fine
     names(v) <- c("A", "B", "C")
     ## checkEquals(v, w)            ## fails because v and w have different names
     checkEqualsNumeric(v, w)        ## passes fine because names are ignored
     x <- rep(1:12, 2)
     y <- rep(0:1, 12)
     res <- list(a=1:3, b=letters, LM=lm(y ~ x))
     res2 <- list(a=seq(1,3,by=1), b=letters, LM=lm(y ~ x))
     checkEquals( res, res2)        ## passes fine
     checkIdentical( res, res)
     checkIdentical( res2, res2)
     ## checkIdentical( res, res2)  ## fails because element 'a' differs in type
     fun <- function(x) {
       if(x)
       {
        stop("stop conditions signaled")
       }
       return()
     }
     checkException(fun(TRUE))      ## passes fine
     ## checkException(fun(FALSE))  ## failure, because fun raises no error
     checkException(fun(TRUE), silent=TRUE)
     ##  special constants
     ##  same behaviour as for underlying base functions
     checkEquals(NA, NA)
     checkEquals(NaN, NaN)
     checkEquals(Inf, Inf)
     checkIdentical(NA, NA)
     checkIdentical(NaN, NaN)
     checkIdentical(-Inf, -Inf)
}


# Add a bug with automatic interval/bin detection on binary variable (all vals get placed in one bin)
test.bin01bug <- function() {

}


get.testDat <- function(nsamp = 100000) {
  `%+%` <- function(a, b) paste0(a, b)
  library(simcausal)
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
  D <- set.DAG(D)
  datO <- sim(D, n = nsamp, rndseed = 12345)
}

get.testDatNet <- function(datO) {
  Kmax <- 1
  # nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"))
  nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"), nFnode = "nF")
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA")
  netind_cl <- NetIndClass$new(nobs = nrow(datO))
  # Define datNetObs:
  datnetW <- DatNet$new(netind_cl = netind_cl, nodes = nodes, VarNodes = nodes$Wnodes, addnFnode = TRUE)$make.sVar(Odata = datO, sVar.object = def_sW)
  datnetA <- DatNet$new(netind_cl = netind_cl, nodes = nodes, VarNodes = nodes$Anode)$make.sVar(Odata = datO, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
  return(list(datNetObs = datNetObs, netind_cl = netind_cl, def_sA = def_sA, def_sW = def_sW, nodes = nodes))
}


test.RegressionClass <- function() {
  # Tests for RegressionClass:
  reg_test1 <- RegressionClass$new(outvar.class = c(gvars$sVartypes$bin, gvars$sVartypes$bin),
                                  outvar = c("A1", "A2"),
                                  predvars = c("W1", "W2"))
  class(reg_test1)
  reg_test1$subset
  model1 <- SummariesModel$new(reg = reg_test1)
  # [1] "Init BinOutModel:"
  # [1] "P(A1|W1,W2)"
  # [1] "Init BinOutModel:"
  # [1] "P(A2|A1,W1,W2)"
  class(model1$getPsAsW.models()$`P(sA|sW).1`$reg)
  ls(model1$getPsAsW.models()$`P(sA|sW).1`$reg)

  reg_test2 <- RegressionClass$new(outvar.class = c(gvars$sVartypes$bin, gvars$sVartypes$bin),
                                  outvar = c("A1", "A2"),
                                  predvars = c("W1", "W2"),
                                  subset = list(quote(A1 == 0)))
  class(reg_test2)
  reg_test2$subset
  model2 <- SummariesModel$new(reg = reg_test2)

  reg_test3 <- RegressionClass$new(outvar.class = c(gvars$sVartypes$cont, gvars$sVartypes$cont),
                                  outvar = c("sA"),
                                  predvars = c("W1", "W2", "W3"),
                                  subset = list(quote(sA==0)))
  reg_test_new <- reg_test3$clone()
  all.equal(reg_test3, reg_test_new)
  class(reg_test3)
  nsamp <- 10000
  datO <- get.testDat(nsamp)
  head(datO)
  nodeobjs <- get.testDatNet(datO)
  model3 <- SummariesModel$new(reg = reg_test3, O.datnetA = nodeobjs$datNetObs$datnetA)
}

test.PoolContRegression <- function() {
  library(data.table)

  reg_test <- RegressionClass$new(outvar.class = c(gvars$sVartypes$cont, gvars$sVartypes$cont),
                                  outvar = c("sA"),
                                  predvars = c("W1", "W2", "W3"),
                                  subset = list(quote(TRUE)))
  # datO <- get.testDat(nsamp = 50000)
  # datO <- get.testDat(nsamp = 100000)
  datO <- get.testDat(nsamp = 10000)
  nodeobjs <- get.testDatNet(datO)
  datNetObs <- nodeobjs$datNetObs
  class(datNetObs) # [1] "DatNet.sWsA" "DatNet"      "R6"

  model3 <- SummariesModel$new(reg = reg_test, O.datnetA = nodeobjs$datNetObs$datnetA)
  # Matrix of all summary measures: (sW,sA)
  head(nodeobjs$datNetObs$mat.sVar); class(nodeobjs$datNetObs$mat.sVar)
  binfit_time <- system.time(
    model3$fit(data = nodeobjs$datNetObs)
  )
  binfit_time
  # for 50K obs:
  # user  system elapsed 
  # 5.136   3.622   8.734 
  # for 10K obs:
  #  user  system elapsed 
  # 0.199   0.049   0.243   

  # [1] "fit (10K)"
  # $coef
  #  Intercept     bin_ID         W1         W2         W3 
  # -2.7756215  0.1553186 -1.0014477 -0.5720651 -0.3339728 
  # [1] "res_DT: "
  #            ID ProbAeqa_long
  #      1:     1    0.97396036
  #      2:     1    0.96971742
  #      3:     1    0.96480811
  #      4:     1    0.95913645
  #      5:     1    0.95259565
  #     ---                    
  # 104496:  9998    0.07668105
  # 104497:  9999    0.93215687
  # 104498:  9999    0.92165035
  # 104499:  9999    0.09032560
  # 104500: 10000    0.06784313
  # [1] "res_DT_short: "
  #           ID    cumprob
  #     1:     1 0.06099655
  #     2:     2 0.06145986
  #     3:     3 0.03836225
  #     4:     4 0.05821479
  #     5:     5 0.07303417
  #    ---                 
  #  9996:  9996 0.05119563
  #  9997:  9997 0.05896735
  #  9998:  9998 0.06414013
  #  9999:  9999 0.07760077
  # 10000: 10000 0.06784313
  # [1] "head(ProbAeqa, 50)"
  #  [1] 0.060996548 0.061459862 0.038362248 0.058214786 0.073034166 0.064140127 0.060658764 0.050023002 0.026039639 0.075325033 0.029168620
  # [12] 0.054538219 0.054031618 0.083549453 0.008653412 0.029594466 0.077600772 0.081220201 0.068319822 0.061459862 0.071357407 0.039453938
  # [23] 0.075325033 0.039007914 0.057871503 0.077600772 0.057871503 0.058967354 0.064140127 0.043973691 0.046655735 0.079794387 0.074434114
  # [34] 0.058967354 0.067843133 0.063492979 0.033237556 0.064138704 0.056974041 0.065426910 0.037236039 0.029168620 0.056974041 0.047226347
  # [45] 0.043973691 0.084256432 0.060173071 0.073034166 0.029168620 0.060183301

}



# NOT IMPLEMENTED
# Test with bivariate normal (as sA)
test.iid.bivNorm.tmlefit <- function() {

  rbivNorm <- function(n, whichbiv, norms, mu, var1 = 1, var2 = 1, rho = 0.7) {
    whichbiv <- whichbiv[1]; var1 <- var1[1]; var2 <- var2[1]; rho <- rho[1]
    sigma <- matrix(c(var1, rho, rho, var2), nrow = 2)
    Scol <- chol(sigma)[, whichbiv]
    bivX <- (Scol[1] * norms[, 1] + Scol[2] * norms[, 2]) + mu
    bivX
  }

  D <- DAG.empty()
  D <-
  D + 
      # W1, W2, W3, A1_0, A2_0 will play the role of sW:
      node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("A1.norm1", distr = "rnorm", mean = 0, sd = 1) +
      node("A1.norm2", distr = "rnorm", mean = 0, sd = 1) +
      node("alpha", t = 0:1, distr = "rconst",
                const = {if(t == 0) {log(0.6)} else {log(1.0)}}) +

      # A1_1, A2_1, Y1_1 will play the role of sA:
      node("A1", t = 0:1, distr = "rbivNorm", whichbiv = t + 1,
                norms = c(A1.norm1, A1.norm2),
                mu = {if (t == 0) {0} else {-0.30 * A1[t-1]}}) +
      node("A2", t = 0:1, distr = "rbern",
                prob = plogis(alpha[t] + log(5)*A1[t] + {if(t == 0) {0} else {log(5)*A2[t-1]}})) +
      node("Y1", t = 1, distr = "rnorm",
                mean = (0.98 * W1 + 0.58 * W2 + 0.33 * W3 +
                        0.98 * A1[t] - 0.37 * A2[t]),
                sd = 1) +

      # outcome:
      node("Y2", t = 1, distr = "rnorm",
                mean = (Y1[t] + 0.98 * W1 + 0.58 * W2 + 0.33 * W3 +
                        0.98 * A1[t] - 0.37 * A2[t]),
                sd = 1)

  D <- set.DAG(D)
  datO <- sim(D, n = nsamp, rndseed = 12345)
  head(datO)

  hist(datO$A1_0)
  hist(datO$A1_1)
  plot(density(datO$A1_0))
  lines(density(datO$A1_1))

  hist(datO$Y1_1)
  hist(datO$Y2_1)
  plot(density(datO$Y1_1))
  lines(density(datO$Y2_1))

  kmax <- 1
  def_sW <- def.sW(W1="W1", W2="W2", W3="W3", A1_0="A1_0", A2_0="A2_0")
  def_sA <- def.sA(A1_0="A1_1", A2_0="A2_1", Y1_1="Y1_1")
  tmlenet_res <- tmlenet(data = datO, Anode = "A2_1", Wnodes = c("W1", "W2", "W3", "A1_0", "A2_0"), Ynode = "Y2_1",
                          Kmax = kmax, IDnode = "ID", NETIDnode = "ID",
                          f_gstar1 = f.A_0,
                          sW = def_sW, sA = def_sA,
                          # Qform = Qform, hform = hform, #gform = gform,  # remove
                          # new way to specify regressions:
                          Qform = "Y2_1 ~ W1 + W2 + W3 + A1_0 + A2_0 + A1_1 + A2_1 + Y1_1",
                          hform = "A1_1 + A2_1 + Y1_1 ~ W1 + W2 + W3 + A1_0 + A2_0",
                          hform.gstar = "A1_1 + A2_1 + Y1_1 ~ W1 + W2 + W3 + A1_0 + A2_0",
                          gform = "A2 ~ W1 + W2 + W3 + A1_0 + A2_0",
                          optPars = list(
                            onlyTMLE_B = TRUE,  # remove
                            # f_g0 = f.A, # tested, works
                            n_MCsims = 10
                          )
                          )

}




## ---------------------------------------------------------------------
# TEST SET FOR DETECTING VECTOR TYPES
## ---------------------------------------------------------------------
test_mat <- as.matrix(data.frame(a = c(0,1,0,0,1), b = rep(5,5), c = c(1,2,3,4,5), d = rnorm(5)))
correct.types <- list(a = gvars$sVartypes$bin, b = gvars$sVartypes$bin, c = gvars$sVartypes$cat, d = gvars$sVartypes$cont)
out.types <- detect.col.types(test_mat)
all.equal(correct.types, out.types)

## ---------------------------------------------------------------------
# TEST SET FOR NetIndClass class
## ---------------------------------------------------------------------
testdf <- data.frame(a = rnorm(100), b = rnorm(100))
nettest <- NetIndClass$new(Odata = testdf, Kmax = 5)
nettest$getNetInd
nettest$getNetInd <- NULL

## ---------------------------------------------------------------------
# TESTING ContinOutModel class and NewSummaryModel.contin constructor
## ---------------------------------------------------------------------
# TEST 1: Binary outvar (treated as continuous). Testing that results match with binary class prediction.
# self <- list()
# self$nbins <- 10L
# self$cats <- ...
# self$sA_nms <- ...
# (binvar <- rbinom(n=100, size = 1, prob = 0.3))
# (int_bylen <- c(0, 0.1, 0.9, 1))
# (int_bymass <- quantile(binvar, prob = c(0, 0.1, 0.9, 1)))
# (ord1 <- findInterval(x = binvar, vec = int_bylen, rightmost.closed = TRUE))
# (ord2 <- findInterval(x = binvar, vec = int_bymass, rightmost.closed = TRUE))
# make.bins_mtx_1(x.ordinal = ord.sAj, self = self)
# obj <- NewSummaryModel.contin()
# obj$fit()

# sAnorm <- normalize(x = rnorm(100000))
# # sAnorm <- normalize(x = datatest$sA)
# (ints_list <- define.intervals(x = sAnorm, nbins = self$nbins))
# xcat_bylen = discretize(x = sAnorm, intervals = ints_list$intbylen) # NOTE: x as categorical for intervals def by equal len 0-1
# xcat_bymass = discretize(x = sAnorm, intervals = ints_list$intbymass) # NOTE x as categorical for intervals def by equal mass 0-1
# # data.frame(sAnorm, xcat_bylen, xcat_bymass)
# hist(xcat_bylen)
# hist(xcat_bymass)
# sort(unique(xcat_bylen))
# sort(unique(xcat_bymass))

# (t1 <- system.time(bins_mat1 <- make.bins_mtx_1(x.ordinal = xcat_bylen)))
# outdf1 <- data.frame(sAnorm, xcat_bylen, bins_mat1)

# (t2 <- system.time(bins_mat2 <- make.bins_mtx_1(x.ordinal = xcat_bymass)))
# outdf2 <- data.frame(sAnorm, xcat_bymass, bins_mat2)

# print(head(outdf1))
# print(head(outdf2))
# head(cbind(bins_mat1, bins_mat2))

# TEST 2: Continuous outvar.
# datatest <- data.frame(
#                       W = rbinom(100, 1, prob=.5), 
#                       A = rbinom(100, 1, prob=.2),
#                       sA = rnorm(100),
#                       Y = rbinom(100, 1, prob=.7), 
#                       nFriends = rbinom(100, 5, prob=0.2)
#                       )
# nodes <- list(Anode = "A", Wnode = "W", Ynode = "Y", nFnode = "nFriends")


## ---------------------------------------------------------------------
# TESTING classes not throw exceptions for fitting empty X_mat design mat and empty Y outcome vector
## ---------------------------------------------------------------------
# testmat <- matrix(rnorm(50), nrow=10, ncol = 5)
# colnames(testmat) <- "A"%+%(1:5)
# Yvals <- rep(1, 10)
# fit <- speedglm.wfit(y = Yvals, X = testmat, family = binomial())$coef
# class(fit)
# attributes(fit)

# testmat_emp <- testmat[FALSE, ]
# nrow(testmat_emp)==0L
# Yvals_emp <- Yvals[FALSE]
# cbind(Yvals_emp, testmat_emp)
# # both throw exceptions for empty design mat:
# speedglm.wfit(y = Yvals_emp, X = testmat_emp, family = binomial())$coef
# glm.fit(x = testmat_emp, y = Yvals_emp, family = binomial())
# logisfit.glmS3 = function(datsum_obj) { # # S3 method for glm binomial family fit, takes DatBinarySummary data object 
# logisfit.speedglmS3 = function(datsum_obj) { # S3 method for speedglm binomial family fit, takes DatBinarySummary data object 





test.continous.sA <- function() {
    # I) Build network vectors: (W, W_netF_1, ..., W_netF_k) for each W in Wnodes by PRE-ALLOCATING netW_full:
    datnetW <- DatNet$new(Odata = data, NetInd_k = NetInd_k, Kmax = k, nodes = node_l, VarNodes = node_l$Wnodes, addnFnode = TRUE)
    # II) APPLY THE SUMMARY MEASURE FUNCTIONS / EXPRESSION TO netW_full to OBTAIN sW columns SELECT ONLY sW columns in hform_g0 and hfrom_gstar or use all? 
    obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms, norm.c.sVars = TRUE)$dat.sVar
    datnetW$def_cbin_intrvls()
    print("Detected intervals: "); print(datnetW$cbin_intrvls)
    print("Detected nbins: "); print(datnetW$nbins)

    oldopts <- tmlenet_options(maxncats = 5, nbins = 10)

    print("No normalization. Binning by mass")
      obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms)$dat.sVar
      print("head(obsdat.sW)"); print(head(obsdat.sW))
      print("Var types: "); str(datnetW$type.sVar)
      print("Contin covars names: "); print(datnetW$names.c.sVar)
      defints1 <- datnetW$def_cbin_intrvls()$cbin_intrvls
      print("No normalization bin intervals by mass: "); print(defints1)
      print("nbins: "); print(datnetW$nbins); 
      print("Testing ordinals with ncats < nbins get nbins = ncats:"); print(datnetW$nbins["nFriends"] < tlmenet:::getopt("nbins"))

    print("No normalization. Binning by equal length")
      intlen_ints <- datnetW$def_cbin_intrvls(bin_bymass = FALSE)$cbin_intrvls
      print("No normalization bin intervals by length: "); print(intlen_ints)
      print("nbins: "); print(datnetW$nbins)

    print("Testing ordinals with ncats > nbins get collapsed into fewer cats:")
    tmlenet_options(nbins = 4)
      obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms)$dat.sVar
      defints2 <- datnetW$def_cbin_intrvls()$cbin_intrvls
      print("New bins with collapsed ordinals: "); print(defints2)
      print("nbins: "); print(datnetW$nbins)
    tmlenet_options(nbins = 10)
    print("Testing normalization:")
      obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms, norm.c.sVars = TRUE)$dat.sVar
      print("head(obsdat.sW)"); print(head(obsdat.sW))
      savedtypes <- datnetW$type.sVar
      print("Var types: "); str(datnetW$type.sVar)
      print("Contin covars names: "); print(datnetW$names.c.sVar)
      savedc.ints <- datnetW$def_cbin_intrvls()$cbin_intrvls
      print("Detected bin intervals: "); print(savedc.ints)
      print("nbins: "); print(datnetW$nbins)

    print("testing overwriting of bin intervals 1:")
      newint <- seq(from = 0, to = 1, length.out = (tlmenet:::getopt("nbins")+1))
      datnetW$def_cbin_intrvls(cbin_intrvls = newint)
      print("Assigned bin intervals 1: "); print(datnetW$cbin_intrvls)
      print("nbins: "); print(datnetW$nbins)
      print("Correctly assigned bin intervals 1: "); print(all.equal(as.vector(datnetW$cbin_intrvls[[1]]), as.vector(datnetW$cbin_intrvls[[2]])))

    print("testing overwriting of bin intervals 2:")
      int1 <-  seq(from = 0, to = 0.5, length.out = (tlmenet:::getopt("nbins")+1))
      int2 <-  seq(from = 0.5, to = 1, length.out = (tlmenet:::getopt("nbins")+1))
      newcbins <- list(nFriends = int1, netW3_sum = int2)
      datnetW$def_cbin_intrvls(cbin_intrvls = newcbins)
      print("Assigned bin intervals 2: "); print(datnetW$cbin_intrvls);
      print("nbins: "); print(datnetW$nbins)
      print("Correctly assigned bin intervals 2: "); print(all.equal(datnetW$cbin_intrvls, newcbins))

    print("testing passing custom type.sVar list:")
      datnetW$make_sVar(names.sVar = sW_nms, type.sVar = gvars$sVartypes$bin, norm.c.sVars = TRUE)
      print("assigned all the same bin types: ");
      print(all(unlist(lapply(datnetW$type.sVar, function(x) x%in%gvars$sVartypes$bin))))
      str(datnetW$type.sVar)
      print("Cont. var names correct: "); print(length(datnetW$names.c.sVar) == 0L)
      print("Redefined c bin intervals correct: "); print(length(datnetW$def_cbin_intrvls()$cbin_intrvls) == 0L)
      print("nbins: "); print(datnetW$nbins)

      datnetW$make_sVar(names.sVar = sW_nms, type.sVar = savedtypes, norm.c.sVars = TRUE)
      print("Re-assigned old types: "); str(datnetW$type.sVar)
      newc.ints <- datnetW$def_cbin_intrvls()$cbin_intrvls
      print("Redefined c bin intervals: "); print(newc.ints)
      print("nbins: "); print(datnetW$nbins)
      print("Redefined c bin ints match old ones: "); print(all.equal(savedc.ints, newc.ints))

    do.call(tmlenet_options, oldopts)
}



notest.evalsubset2dfs <- function() {
  # Testing eval'ing logical expressions in envirs of two data.frames with different number of rows.
  mat1 <- matrix(c(0,1,1,1,0), nrow = 5, ncol = 2)
  colnames(mat1) <- c("A1", "A2")
  mat2 <- matrix(c(2,3), nrow = 15, ncol = 4)
  colnames(mat2) <- c("B1", "B2", "B3", "B4")
  # using env as lists from df(mat1) and df(mat2) put together (mat1 and mat2 have diff no. of rows)
  (testenv1 <- c(data.frame(mat1), data.frame(mat2)))
  # using env as one data.frame put together (result has the same no. of rows)
  (testenv2 <- cbind(data.frame(mat1), data.frame(mat2)))
  subsetexpr1 <- parse(text = "A1 == 1")[[1]]
  eval(subsetexpr1, envir = testenv1, enclos = baseenv())
  # [1] FALSE  TRUE  TRUE  TRUE FALSE 
  subsetexpr2 <- parse(text = "B1 == 2")[[1]]
  eval(subsetexpr2, envir = testenv1, enclos = baseenv())
  # [1]  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE
  subsetexpr3 <- parse(text = "(B1 == 2) & (A1 == 1)")[[1]] # automatically repeates the result of the first logical expr (B1 == 2)
  (resenv1 <- eval(subsetexpr3, envir = testenv1, enclos = baseenv()))  
  # [1] FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE
  (resenv2 <- eval(subsetexpr3, envir = testenv2, enclos = baseenv()))
  # [1] FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE
  all.equal(resenv1, resenv2) # works as expected
  # [1] TRUE 
}

test.NetIndClass <- function() {
	# ----------------------------------------------------------------------------------------
	# TESTING NetIndClass CLASS
	# ----------------------------------------------------------------------------------------
	k <- 2
	dftestW <- data.frame(W = as.integer(c(6,7,8,9,10))) # W_netF1 = rep(6,5), W_netF2 = rep(8,5)
	dftestA <- data.frame(A = as.integer(c(1,2,3,4,5))) # , A_netF1 = rep(1,5), A_netF2 = rep(3,5)
	class(dftestW$W)
	class(dftestA$A)

	NET_id <- c(rep("1 3", nrow(dftestW)-1), "1")
	class(dftestW$A)
	dftest1 <- data.frame(dftestW, dftestA, NETID = NET_id)
	is.factor(dftest1$NETID)
	netindcl <- NetIndClass$new(Odata = dftest1, Kmax = k, NETIDnode = "NETID")

	netindcl$NetInd_k
	"matrix" %in% class(netindcl$NetInd_k)
	"integer" %in% class(netindcl$NetInd_k[,1])

	dftest2 <- data.frame(dftestW, dftestA, NETID = NET_id, stringsAsFactors = FALSE)
	is.character(dftest2$NETID)
	netindcl <- NetIndClass$new(Odata = dftest2, Kmax = k, NETIDnode = "NETID")

	netindcl$NetInd_k
	"matrix" %in% class(netindcl$NetInd_k)
	"integer" %in% class(netindcl$NetInd_k[,1])
}



# TESTING sVar expressions parser:
test.Define_sVar <- function() {
  # ----------------------------------------------------------------------------------------
  # TEST DATA:
  # ----------------------------------------------------------------------------------------
  `%+%` <- function(a, b) paste0(a, b)
  k <- 2
  dftestW <- data.frame(W = as.integer(c(6,7,8,9,10)))
  dftestA <- data.frame(A = as.integer(c(1,2,3,4,5))) 
  NET_id <- c(rep("1 3", nrow(dftestW)-1), "1 ")
  # NET_id <- c(rep("1 3", nrow(dftestW)-1), "1")
  class(dftestW$W)
  class(dftestA$A)

  dfnet <- data.frame(dftestW, W_netF1 = rep(6,5), W_netF2 = c(rep(8,4), NA), dftestA, A_netF1 = rep(1,5), A_netF2 = c(rep(3,4), NA))
  dfnet

  # ----------------------------------------------------------------------------------------
  # TESTING sVar expressions parser
  # ----------------------------------------------------------------------------------------
  dftest <- data.frame(dftestW, dftestA, NETID = NET_id)
  NetInd_cl <- NetIndClass$new(Odata = dftest, Kmax = k, NETIDnode = "NETID")
  NetInd_cl$Kmax
  NetInd_cl$NetInd_k
  NetInd_cl$nF
  NetInd_cl$mat.nF

  # **** Example TESTING Kmax substitute ****
  defsVar.expr0 <- def.sW.g0(sA.1 = A[[Kmax]])
  (evaled.sVar.expr0 <- defsVar.expr0$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  (evaled.sVar.expr0 <- defsVar.expr0$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl, addnFnode = TRUE)$mat.sVar)

  # Example 0.
  defsVar.expr0 <- def.sW.g0(sA.1 = A)
  all(as.vector(evaled.sVar.expr0) == dftest$A)

  (evaled.sVar.expr0 <- defsVar.expr0$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  defsVar.expr0 <- def.sW.g0(A)
  (evaled.sVar.expr0 <- defsVar.expr0$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  defsVar.expr0 <- def.sW.g0(A[[0]])
  (evaled.sVar.expr0 <- defsVar.expr0$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)  

  defsVar.expr0 <- def.sW.g0(A[[0:Kmax]])
  (evaled.sVar.expr0 <- defsVar.expr0$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)  

  class(defsVar.expr0)
  class(evaled.sVar.expr0)
  is.matrix(evaled.sVar.expr0)

  # Example 1.
  defsVar.expr1 <- def.sW.g0(sA.1 = rowSums(A[[0:k]]))
  (evaled.sVar.expr1 <- defsVar.expr1$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  # w/ NA for missing vars:
  (evaled.sVar.expr1 <- defsVar.expr1$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl, misXreplace = gvars$misval)$mat.sVar)
  is.na(evaled.sVar.expr1[5,1])
  (evaled.sVar.expr1 <- defsVar.expr1$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl, misXreplace = 999)$mat.sVar)
  evaled.sVar.expr1[5,1]==1005

  # Example 2. Using a variable to pass sVar expression.
  (testexpr_call <- quote(rowSums(A[[0:k]])))
  # (defsVar.expr2 <- def.sW.g0(W = testexpr_call)) # doesn't work
  defsVar.expr2 <- def.sW.g0(sA.1 = eval(testexpr_call))
  class(defsVar.expr2$sVar.exprs[[1]])
  (evaled.sVar.expr2 <- defsVar.expr2$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  res1 <- as.integer(c(5, 6, 7, 8, 6))

  evaled.sVar.expr1 <- defsVar.expr1$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar
  all.equal(evaled.sVar.expr1, evaled.sVar.expr2)
  all(res1 == as.vector(evaled.sVar.expr1))

  # Example 3. Generate a matrix of sVar[1], ..., sVar[j] from one sVar expression.
  defsVar.expr1 <- def.sW.g0(W = W[[0:k]])
  (evaled.sVar.expr1 <- defsVar.expr1$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)

  defsVar.expr1 <- def.sW.g0(W[[0:k]])
  (evaled.sVar.expr1 <- defsVar.expr1$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)

  is.matrix(evaled.sVar.expr1)
  class(evaled.sVar.expr1)  # [1] "matrix"
  dim(evaled.sVar.expr1)  # [1] 5 3
  all(evaled.sVar.expr1[,1] == dftest$W)
  defsVar.expr2 <- def.sW.gstar(W = W[[0:k]])
  class(defsVar.expr2)
  (evaled.sVar.expr2 <- defsVar.expr2$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  all.equal(evaled.sVar.expr1, evaled.sVar.expr2)

  # Example 4a. Generate a matrix of sVar[1], ..., sVar[j] from one sVar expression that is a combination of different Vars in Odata.
  defsVar.expr <- def.sA(sA.1 = W[[0:k]] + rowSums(A[[1:k]]))
  (evaled.sVar.expr <- defsVar.expr$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  class(evaled.sVar.expr)
  colnames(evaled.sVar.expr)

  testres1_cl <- def.sA(netW = W[[0:k]])
  (evaled.testres1 <- testres1_cl$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  testres2_cl <- def.sA(sA.1 = rowSums(A[[1:k]]))
  evaled.testres2 <- testres2_cl$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar
  all((evaled.testres1 + as.vector(evaled.testres2)) == evaled.sVar.expr)

  # Example 4b. Generate a matrix of sVar[1], ..., sVar[j] from one sVar expression that is a combination of different Vars in Odata.
  defsVar.expr <- def.sW.g0(W = "W[[0:k]] + rowSums(A[[1:k]])")
  class(defsVar.expr$sVar.exprs[["W"]])
  defsVar.expr$sVar.exprs[["W"]]
  (evaled.sVar.expr2 <- defsVar.expr$parse.sVar(data.df = dftest,  NetInd_cl = NetInd_cl)$mat.sVar)
  all(evaled.sVar.expr2==evaled.sVar.expr)

  # Example 5. sum of prod of netA and netW:
  defsVar.expr <- def.sA(sumAWnets = rowSums(A[[1:k]] * W[[1:k]]))
  (evaled.sVar.expr <- defsVar.expr$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
  all(as.integer(as.vector(evaled.sVar.expr)) == c(30,30,30,30,6))

  # Example 6. More than one summary measure
  defsVar.expr <- def.sA(A = A, sumAnets = rowSums(A[[1:k]]), sumAWnets = rowSums(A[[1:k]] * W[[1:k]]))
  (evaled.sVar.expr <- defsVar.expr$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)

  # Example 7. No names
  defsVar.expr <- def.sA(A, sumAnets = rowSums(A[[1:k]]), sumAWnets = rowSums(A[[1:k]] * W[[1:k]]))
  (evaled.sVar.expr <- defsVar.expr$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)

  defsVar.expr <- def.sA(A[[0:Kmax]], sumAnets = rowSums(A[[1:k]]), sumAWnets = rowSums(A[[1:k]] * W[[1:k]]))
  (evaled.sVar.expr <- defsVar.expr$parse.sVar(data.df = dftest, NetInd_cl = NetInd_cl)$mat.sVar)
}




test.bugfixes <- function() {

}





test.plotting <- function() {

}





test.node <- function() {

}