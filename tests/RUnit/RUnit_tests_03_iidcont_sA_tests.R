# ---------------------------------------------------------------------------------
# TEST SET 2. TESTS FOR FITTING CONTINUOUS EXPOSURE sA IN IID DATA
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by  binning, conditional on covariates
# Overall exposure g0 (sA) is a mixture of 3 normals,
# individual exposure is normal with mu for each observation being a function of (W1,W2,W3), sd = 1;
# ---------------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA <- function(x) all(is.na(x))

# ---------------------------------------------------------------------------------
# Test 1. Directly fit a joint density for sA, sA2 ~ W, for sA - continuous (with speedglm and glm.fit)
# ---------------------------------------------------------------------------------
  ## helper fun
  get.density.sAdat <- function(nsamp = 100000) {
    require(simcausal)
    D <- DAG.empty()
    D <-
    D + node("W1", distr = "rbern", prob = 0.5) +
        node("W2", distr = "rbern", prob = 0.3) +
        node("W3", distr = "rbern", prob = 0.3) +
        node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
        node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
    D <- set.DAG(D, n.test = 10)
    datO <- sim(D, n = nsamp, rndseed = 12345)
  }
  ## helper fun
  get.setW.sAdat <- function(setWvals, nsamp = 100000) {
    require(simcausal)
    W1val <- setWvals[1]; W2val <- setWvals[2]; W3val <- setWvals[3]
    D <- DAG.empty()
    D <-
    D + node("W1", distr = "rconst", const = .(W1val)) +
        node("W2", distr = "rconst", const = .(W2val)) +
        node("W3", distr = "rconst", const = .(W3val)) +
        node("sA", distr = "rnorm", mean = (0.98 * W1 + 0.58 * W2 + 0.33 * W3), sd = 1)
    D <- set.DAG(D, n.test = 10)
    datWset <- sim(D, n = nsamp, rndseed = 12345)
    setWmat <- as.matrix(data.frame(W1 = W1val, W2 = W2val, W3 = W3val, sA = seq(-4, 4, by = 0.2)))
    return(list(setWsA = datWset$sA, setWmat = setWmat))
  }
  ## helper fun
  def.nodeojb <- function(datO) {
    Kmax <- 1
    nodes <- list(Anodes = "sA", Wnodes = c("W1", "W2", "W3"), nFnode = "nF")
    sW <- def_sW(W1 = "W1", W2 = "W2", W3 = "W3")
    sA <- def_sA(sA = "sA")
    netind_cl <- simcausal::NetIndClass$new(nobs = nrow(datO))
    # Define datNetObs:
    OdataDT_R6 <- OdataDT$new(Odata = datO, nFnode = "nF", iid_data_flag = FALSE)
    datnetW <- DatNet$new(Odata = OdataDT_R6, netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = OdataDT_R6, sVar.object = sW)
    checkTrue(tmlenet:::is.DatNet(datnetW))
    datnetA <- DatNet$new(Odata = OdataDT_R6, netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = OdataDT_R6, sVar.object = sA)
    datNetObs <- DatNet.sWsA$new(Odata = OdataDT_R6, datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
    return(list(datNetObs = datNetObs, netind_cl = netind_cl, sA = sA, sW = sW, nodes = nodes))
  }

test.sampleA <- function() {
  nsamp <- 10000
  datO <- get.density.sAdat(nsamp)
  nodeobjs <- def.nodeojb(datO)
  testm.sW <- nodeobjs$sW$eval.nodeforms(data.df = datO, netind_cl = nodeobjs$netind_cl)
  testm.sA <- nodeobjs$sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs$netind_cl)

  DatNet_object <- nodeobjs$datNetObs
  head(DatNet_object$dat.sWsA)
  subset_mat <- DatNet_object$get.dat.sWsA(rowsubset = TRUE, covars = c("W1","W2"))

  # Define est_params_list:
  reg.sVars <- list(outvars = c("sA"), predvars = c("W1", "W2", "W3"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- nodeobjs$datNetObs$datnetA$type.sVar[reg.sVars$outvars]

  # Put all est_params in RegressionClass (regression with speedglm package)
  regclass <- RegressionClass$new(outvar.class = sA_class,
                                  outvar = reg.sVars$outvars,
                                  predvars = reg.sVars$predvars,
                                  subset = subset_vars,
                                  nbins=50)

  summeas.g0 <- SummariesModel$new(reg = regclass, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)
  summeas.g0$predict(newdata = nodeobjs$datNetObs)  # summeas.g0$sA_nms

  ## resample the outcome (sA), conditional on observed predictors, using the current conditional density fit for P(sA|X)
  set.seed(123456)
  resampA_20times <- lapply(1:20, function(i) summeas.g0$sampleA(newdata = nodeobjs$datNetObs))
  resampA_combined <- unlist(resampA_20times)
  resampA_average <- rowMeans(do.call("cbind", resampledA_20times))
  summary(resampledA_average - datO[["sA"]])
  print(mean(resampledA_average) - mean(datO[["sA"]]))
  checkTrue(mean(resampledA_average) - mean(datO[["sA"]]) < 0.005)

  # Get P(sA|W) for the observed data (W,sA):
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  ## ---------------------------------------------------------------------------------------------------------
  ## Plot density of the entire observed outcome vs. resampled outcomes (resampled 20 times and then combined into one long vector)
  ## ---------------------------------------------------------------------------------------------------------
  set.seed(12345)
  dens_A <- density(datO[["sA"]])
  dens_resampA <- density(resampA_combined)

  ## plot the density fit with observed A vs. resampled A:
  plot(dens_A)
  lines(dens_resampA, col = "blue")

  ## numerically compare the two density fits:
  common_x <- intersect(round(dens_A[["x"]], 2), round(dens_resampA[["x"]], 2))
  idx_A <- which(round(dens_A[["x"]], 2) %in% common_x)
  idx_resampA <- which(round(dens_resampA[["x"]], 2) %in% common_x)
  length(idx_resampA) == length(idx_A)
  ## MSE on the density fit of the resampled A values (compared to density fit of the observed A values)
  MSE_resampA <- mean((dens_A[["y"]][idx_A] - dens_resampA[["y"]][idx_resampA])^2)
  print("MSE for observed density fit vs. resampled: " %+% MSE_resampA)
  checkTrue(MSE_resampA < 1e-5)

  ## ---------------------------------------------------------------------------------------------------------
  ## Same as above but conditional on some FIXED values of predictors (W1,W2,W3)
  ## + adding the predicted discretized probs conditional on same W's
  ## ---------------------------------------------------------------------------------------------------------
  setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
  subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])
  setWdat_res <- get.setW.sAdat(setWvals, nsamp)

  ## plot densitity first:
  plot(density(setWdat_res$setWsA))
  lines(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")
  set.seed(12345)
  resampA <- summeas.g0$sampleA(newdata = nodeobjs$datNetObs)
  ## add the density of the resamp outcomes:
  lines(density(resampA[subs]), col = "blue")

}

test.simple.fit.density.sA <- function() {
  nsamp <- 10000
  datO <- get.density.sAdat(nsamp)
  nodeobjs <- def.nodeojb(datO)
  testm.sW <- nodeobjs$sW$eval.nodeforms(data.df = datO, netind_cl = nodeobjs$netind_cl)
  testm.sA <- nodeobjs$sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs$netind_cl)

  DatNet_object <- nodeobjs$datNetObs
  head(DatNet_object$dat.sWsA)
  subset_mat <- DatNet_object$get.dat.sWsA(rowsubset = TRUE, covars = c("W1","W2"))
  # check error is thrown when covariate doesn't exist:
  checkException(subset_mat_error <- DatNet_object$get.dat.sWsA(rowsubset = TRUE, covars = c("W1","W4", "W5")))

  # Define est_params_list:
  reg.sVars <- list(outvars = c("sA"), predvars = c("W1", "W2", "W3"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- nodeobjs$datNetObs$datnetA$type.sVar[reg.sVars$outvars]

  # Put all est_params in RegressionClass (regression with speedglm package)
  print("fitting h_gN based equal.len intervals (default) and speedglm (default): ")
  regclass <- RegressionClass$new(outvar.class = sA_class,
                                  outvar = reg.sVars$outvars,
                                  predvars = reg.sVars$predvars,
                                  subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)

  # Test the coef and summary functions for binoutmodel class:
  out_ContinSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(tmlenet:::coef.BinOutModel(out_BinModels[[1]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[2]]))

  # (intrvls <- out_ContinSummaryModel$intrvls)
  # (intrvls.width <- diff(intrvls))
  # length(intrvls.width)
  # ord.sVar <- nodeobjs$datNetObs$discretize.sVar(name.sVar = "sA", intervals = out_ContinSummaryModel$intrvls)
  # ord.sVar_bw <- intrvls.width[ord.sVar]
  # print(head(cbind(sA = nodeobjs$datNetObs$dat.sVar[, "sA"], ord.sVar, bw = ord.sVar_bw, nodeobjs$datNetObs$dat.bin.sVar), 5))
  # print("freq count for transformed ord.sVar: "); print(table(ord.sVar))
  # plot(density(ord.sVar))
  # hist(ord.sVar)
  summeas.g0$predict(newdata = nodeobjs$datNetObs)  # summeas.g0$sA_nms
  # Get P(sA|W) for the observed data (W,sA):
  # SHOULD BE SIMILAR TO THE OBSERVED DENSITY OF s.A (but discretized)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  print("h_gN fit under speedglm: " %+% mean(h_gN)) # [1] 0.2718823
  checkTrue(abs(mean(h_gN)-0.2718823) < 10^-4)
  # ---------------------------------------------------------------------------------------------------------
  # Plot predicted discretized probs conditional on some values of W's
  # ---------------------------------------------------------------------------------------------------------
  setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
  subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])
  sum(subs)
  setWdat_res <- get.setW.sAdat(setWvals, nsamp)

  # plot densitity first:
  plot(density(setWdat_res$setWsA))
  lines(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")

  # plot predicted vals first:
  # plot(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")
  # lines(density(setWdat_res$setWsA))

  resampA <- summeas.g0$sampleA(newdata = nodeobjs$datNetObs)
  lines(density(resampA[subs]), col = "blue")

  # ---------------------------------------------------------------------------------------------------------
  # Plot all predicted discretized probs together (without conditioning on particular subset of W's)
  # ---------------------------------------------------------------------------------------------------------
  # plot(datO[,"sA"], h_gN, type = "p", cex = .3, col = "red")
  # lines(density(datO[,"sA"]))

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but doing regressions with stats::glm.fit ****
  # ---------------------------------------------------------------------------------------------------------
  print("fitting h_gN based equal.len intervals (default) and glm.fit: ")
  regclass.gml <- RegressionClass$new(useglm = TRUE, outvar.class = sA_class,
                                      outvar = reg.sVars$outvars,
                                      predvars = reg.sVars$predvars,
                                      subset = subset_vars)
  summeas.g0.glm <- SummariesModel$new(reg = regclass.gml, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0.glm$fit(data = nodeobjs$datNetObs)

  # Test the coef and summary functions for binoutmodel class:
  out_ContinSummaryModel <- summeas.g0.glm$getPsAsW.models()$`P(sA|sW).1`
  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(tmlenet:::coef.BinOutModel(out_BinModels[[1]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[2]]))

  summeas.g0.glm$predict(newdata = nodeobjs$datNetObs)
  h_gN.glm <- summeas.g0.glm$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  print("h_gN fit under glm: " %+% mean(h_gN.glm)) # [1] 0.2718823
  checkTrue(abs(mean(h_gN.glm)-0.2718823) < 10^-4)
  # plot densitity first:
  # plot(density(setWdat_res$setWsA))
  # lines(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but doing binnning by mass intervals & regressions with speedglm ****
  # ---------------------------------------------------------------------------------------------------------
  # options(tmlenet.verbose = TRUE)
  # gvars$verbose <- TRUE  # set to TRUE after loading all package functions to print all output
  # ls(gvars)

  print("fitting h_gN based on bin_bymass = TRUE and speedglm (default): ")
  regclass.binmass <- RegressionClass$new(useglm = FALSE,
                                          bin_bymass = TRUE,
                                          max_nperbin = 1000,
                                          outvar.class = sA_class,
                                          outvar = reg.sVars$outvars,
                                          predvars = reg.sVars$predvars,
                                          subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass.binmass, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)
  summeas.g0$predict(newdata = nodeobjs$datNetObs)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  mean(h_gN) # [1] 0.2668144
  checkTrue(abs(mean(h_gN)-0.2668144) < 10^-4)

  # get get the observed data:
  nodeobjs$datNetObs$dat.sVar

  # plot densitity first:
  # plot(density(setWdat_res$setWsA))
  # lines(datO[subs][["sA"]], h_gN[subs], type = "p", cex = .3, col = "red")

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but binning using "dhist" function & regressions with speedglm ****
  # ---------------------------------------------------------------------------------------------------------
  print("fitting h_gN based on dhist intervals and speedglm (default): ")
  regclass.bidhist <- RegressionClass$new(useglm = FALSE,
                                          bin_bymass = FALSE,
                                          bin_bydhist = TRUE,
                                          max_nperbin = 1000,
                                          outvar.class = sA_class,
                                          outvar = reg.sVars$outvars,
                                          predvars = reg.sVars$predvars,
                                          subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass.bidhist, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)

  out_ContinSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  print("Intervals for dhist: ")
  print(out_ContinSummaryModel$intrvls)
  # [1] -1003.5240566    -3.5240566    -1.9562682    -0.8766233    -0.2052710     0.2963970     0.7404754     1.2078705
  # [9]     1.6959939     2.3455552     3.3931981     4.9362651  1004.9362651

  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(tmlenet:::coef.BinOutModel(out_BinModels[[1]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[2]]))

  summeas.g0$predict(newdata = nodeobjs$datNetObs)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  print("regclass.bidhist mean(h_gN) under speedglm: " %+% mean(h_gN))
  # [1] 0.276875
  # [1] 0.27687500785635
  checkTrue(abs(mean(h_gN)-0.276875) < 10^-4)
}


# ---------------------------------------------------------------------------------------------------------
# Test 2. Running iid TMLE fit for continous sA
# TMLE for causal effect in i.i.d. data with continuous exposure under continuous stochastic intervention;
# intervention g.star is defined by shifting the normal density of observed sA until g.star/g.0 >= 10,
# then its truncated to be equal to g.0
# ---------------------------------------------------------------------------------------------------------
# Run one TMLE simulation for iid data sampled from get.iid.densityOdat, estimating psi0 under trunced g.star
# ---------------------------------------------------------------------------------------------------------
get.iid.densityOdat <- function(nsamp = 100000, rndseed = NULL, trunc = 10, shift = 2) {
  require(simcausal)
  # nsamp = 100000
  # rndseed = 12345
  trunc <- trunc
  shift <- shift
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("shift",  distr = "rconst", const = .(shift)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("r.obs.sA",  distr = "rconst", const = exp(shift * (sA - sA.mu - shift / 2))) +
      node("trunc.c",  distr = "rconst", const = .(trunc)) +
      node("untrunc.sA.gstar",  distr = "rconst", const = sA + shift) +
      # node("p.gstar.sA", distr = "rconst", const = (1/sqrt(2*.(pi))) * exp((-1/2) * (untrunc.sA.gstar - sA.mu)^2)) +
      # node("p.gstar.sA.gstar", distr = "rconst", const = (1/sqrt(2*.(pi))) * exp((-1/2) * (untrunc.sA.gstar - (sA.mu + shift))^2)) +
      node("r.new.sA",  distr = "rconst", const = exp(shift * (untrunc.sA.gstar - sA.mu - shift / 2))) +
      node("trunc.sA.gstar",  distr = "rconst", const = ifelse(r.new.sA > trunc.c, sA, untrunc.sA.gstar)) +
      node("probY", distr = "rconst", const = plogis(-0.45 * sA - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y", distr = "rbern", prob = plogis(-0.45 * sA - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("probY.gstar", distr = "rconst", const = plogis(-0.45 * trunc.sA.gstar - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y.gstar", distr = "rbern", prob = plogis(-0.45 * trunc.sA.gstar - 0.5 * W1 - 0.58 * W2 - 0.33 * W3))

  D <- set.DAG(D, n.test = 10)
  datO <- sim(D, n = nsamp, rndseed = rndseed)

  head(datO, 100)
  print("mean(datO$Y): " %+% mean(datO$Y));
  # [1] 0.31625
  psi0 <- mean(datO$Y.gstar)
  print("psi0: " %+% psi0)
  # [1] 0.22378

  return(list(psi0 = psi0, datO = datO))
}

# f.gstar,
run.1sim.tmlenet <- function(nsamp, psi0, Qform, newA.gstar, trunc = 10, shift = 2, n_MCsims = 10) {
  datO <- get.iid.densityOdat(nsamp = nsamp, rndseed = NULL, trunc = trunc, shift = shift)$datO
  datO <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  print("head(datO)"); print(head(datO))
  print("summary(datO)"); print(summary(datO))

  Kmax <- 1
  sW <- def_sW(W1 = "W1", W2 = "W2", W3 = "W3")
  sA <- def_sA(sA = "sA")

  # browser()

  tmlenet_res <- tmlenet(data = datO,
                          # Anodes = "sA",
                          Ynode = "Y",
                          Kmax = Kmax,
                          # nFnode = NULL,
                          # f_gstar1 = f.gstar,
                          intervene1.sA = newA.gstar,
                          sW = sW, sA = sA,
                          Qform = Qform,
                          hform.g0 = "sA ~ W1 + W2 + W3",
                          hform.gstar = "sA ~ W1 + W2 + W3",
                          optPars = list(
                            bootstrap.var = FALSE, n.bootstrap = 10
                          ))
                          # correct Q:
                          # Qform = "Y ~ W1 + W2 + W3 + sA",
                          # misspecified Q:
                          # Qform = "Y ~ W3 + sA",

# optPars = list(n_MCsims = n_MCsims)

  CIs <- tmlenet_res$EY_gstar1$IC.CIs
  (tmle_B.CI <- CIs[rownames(CIs)%in%"tmle_B",])
  (h_iptw.CI <- CIs[rownames(CIs)%in%"h_iptw",])
  cover.tmle_B <- ((psi0 <= tmle_B.CI[2]) && (psi0 >= tmle_B.CI[1]))
  cover.h_iptw <- ((psi0 <= h_iptw.CI[2]) && (psi0 >= h_iptw.CI[1]))

  est_mat <- tmlenet_res$EY_gstar1$estimates
  est <- as.vector(est_mat)
  names(est) <- rownames(est_mat)
  cover <- apply(CIs, 1, function(row) ((psi0 <= row[2]) && (psi0 >= row[1])))
  cover2 <- c(tmle_B = cover.tmle_B, h_iptw = cover.h_iptw)
  # print("ests"); print(est)
  # print("CIs"); print(CIs)
  # print("cover"); print(cover)
  # print("cover2"); print(cover2)
  return(list(est = est, cover = cover))
}

# Function that returns a stochastic intervention function intervening on sA, for given shift
# create_f.gstar <- function(shift, trunc) {
#   shift <- shift
#   trunc <- trunc
#   f.gstar <- function(data, ...) {
#     print("shift: " %+% shift)
#     sA.mu <- 0.98 * data[,"W1"] + 0.58 * data[,"W2"] + 0.33 * data[,"W3"]
#     untrunc.sA <- rnorm(n = nrow(data), mean = sA.mu + shift, sd = 1)
#     r.new.sA <- exp(shift * (untrunc.sA - sA.mu - shift / 2))
#     trunc.sA <- ifelse(r.new.sA > trunc, untrunc.sA - shift, untrunc.sA)
#     return(trunc.sA)
#   }
#   return(f.gstar)
# }

# NOTE:ADD THIS TO AN EXAMPLE OF STOCHASTIC INTERVENTION:
test.onesim.iid.tmlefit <- function() {
  # ---------------------------------------------------------------------------------------------------------
  trunc <- 10
  shift <- 1
  nsamp <- 10000
  # nsamp <- 2000
  # ---------------------------------------------------------------------------------------------------------
  # # get true psi.0:
  # datO <- get.iid.densityOdat(nsamp = 100000, rndseed = 12345, trunc = trunc, shift = shift)
  # psi0 <- datO$psi0
  # print("psi0: " %+% psi0)
  # [1] "psi0: 0.239584" (shift=1)
  # [1] "psi0: 0.2241085" (shift=2)

  # ---------------------------------------------------------------------------------------------------------
  # Misspecified Q:
  # ---------------------------------------------------------------------------------------------------------
  set.seed(33556)
  Qform.mis <- "Y ~ W3 + sA" # misspecified Q:
  # f.gstar <- create_f.gstar(shift = shift, trunc = trunc)

  shift <- shift
  trunc <- trunc

  newA.gstar <-  def_new_sA(sA =
    ifelse(exp(shift * (sA + shift - (0.98*W1 + 0.58*W2 + 0.33*W3) - shift/2)) > trunc,
          sA,
          sA + shift))

  res <- run.1sim.tmlenet(nsamp = nsamp, psi0 = 0, Qform = Qform.mis,
                          # f.gstar = f.gstar,
                          newA.gstar = newA.gstar,
                          trunc = trunc, shift = shift
                          # n_MCsims = 10
                          )
  res$est

# [1] "mean(datO$Y): 0.3156"
# [1] "psi0: 0.2419"
# [1] "head(datO)"
#   W1 W2 W3        sA Y
# 1  1  1  0 1.1530363 0
# 2  1  0  0 0.3955352 1
# 3  1  1  0 1.5730009 0
# 4  1  1  0 0.8500365 0
# 5  0  0  1 1.6758910 0
# 6  1  0  0 2.5997302 0
# [1] "summary(datO)"
#        W1               W2              W3               sA                  Y
#  Min.   :0.0000   Min.   :0.000   Min.   :0.0000   Min.   :-3.754226   Min.   :0.0000
#  1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:-0.009311   1st Qu.:0.0000
#  Median :0.0000   Median :0.000   Median :0.0000   Median : 0.750215   Median :0.0000
#  Mean   :0.4976   Mean   :0.306   Mean   :0.3056   Mean   : 0.766478   Mean   :0.3156
#  3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.: 1.547738   3rd Qu.:1.0000
#  Max.   :1.0000   Max.   :1.000   Max.   :1.0000   Max.   : 5.280415   Max.   :1.0000

  #      tmle    h_iptw     gcomp
  # 0.2357986 0.2439521 0.2130597
  # [1] "new MC.ests mat: "
  #         estimate
  # TMLE   0.2357986
  # h_IPTW 0.2439521
  # MLE    0.2130597
  # test 1:
  checkTrue(abs(res$est["tmle"]-0.2357986) < 10^-4)
  # test 2:
  checkTrue(abs(res$est["h_iptw"]-0.2439521) < 10^-4)
  # test 3:
  checkTrue(abs(res$est["gcomp"]-0.2130597) < 10^-4)

  # ---------------------------------------------------------------------------------------------------------
  # Correct Q:
  # ---------------------------------------------------------------------------------------------------------
  set.seed(23466)
  Qform.corr <- "Y ~ W1 + W2 + W3 + sA" # correct Q:
  res2 <- run.1sim.tmlenet(nsamp = nsamp, psi0 = 0, Qform = Qform.corr,
                          # f.gstar = f.gstar,
                          newA.gstar = newA.gstar,
                          trunc = 10, shift = shift,
                          n_MCsims = 10)
  res2$est
  # [1] "new MC.ests mat: "
  #         estimate
  # TMLE   0.2394423
  # h_IPTW 0.2445226
  # MLE    0.2429667
  #      tmle    h_iptw     gcomp
  # 0.2394423 0.2445226 0.2429667
  # test 1:
  checkTrue(abs(res2$est["tmle"]-0.2394423) < 10^-4)
  # test 2:
  checkTrue(abs(res2$est["h_iptw"]-0.2445226) < 10^-4)
  # test 3:
  checkTrue(abs(res2$est["gcomp"]-0.2429667) < 10^-4)

  # ---------------------------------------------------------------------------------------------------------
  # Correct Q w/ glm.fit:
  # ---------------------------------------------------------------------------------------------------------
  old_opts <- tmlenet_options(useglm = FALSE)
  print_tmlenet_opts()
  set.seed(23466)
  Qform.corr <- "Y ~ W1 + W2 + W3 + sA" # correct Q:
  res3 <- run.1sim.tmlenet(nsamp = nsamp, psi0 = 0, Qform = Qform.corr,
                          # f.gstar = f.gstar,
                          newA.gstar = newA.gstar,
                          trunc = 10, shift = shift,
                          n_MCsims = 10)
  res3$est
  # [1] "new MC.ests mat: "
  #         estimate
  # TMLE   0.2394423
  # h_IPTW 0.2445226
  # MLE    0.2429667
  #      tmle    h_iptw     gcomp
  # 0.2394423 0.2445226 0.2429667
  # test 1:
  checkTrue(abs(res2$est["tmle"]-0.2394423) < 10^-4)
  # test 2:
  checkTrue(abs(res2$est["h_iptw"]-0.2445226) < 10^-4)
  # test 3:
  checkTrue(abs(res2$est["gcomp"]-0.2429667) < 10^-4)

}
