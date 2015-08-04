# ---------------------------------------------------------------------------------
# TESTS FOR IID DATA WITH CONTINUOUS EXPOSURE sA
# ---------------------------------------------------------------------------------


`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))

# ---------------------------------------------------------------------------------
# Fitting continuous exposure by  binning, conditional on covariates
# Overall exposure g0 (sA) is a mixture of 3 normals, individual exposure is normal with mu for each observation being a function of (W1,W2,W3), sd = 1;
# ---------------------------------------------------------------------------------
get.setW.sAdat <- function(setWvals, nsamp = 100000) {
  `%+%` <- function(a, b) paste0(a, b)
  W1val = setWvals[1]
  W2val = setWvals[2]
  W3val = setWvals[3]
  library(simcausal)
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rconst", const = .(W1val)) +
      node("W2", distr = "rconst", const = .(W2val)) +
      node("W3", distr = "rconst", const = .(W3val)) +
      node("sA", distr = "rnorm", mean = (0.98 * W1 + 0.58 * W2 + 0.33 * W3), sd = 1)
  D <- set.DAG(D)
  datWset <- sim(D, n = nsamp, rndseed = 12345)
  setWmat <- as.matrix(data.frame(W1 = W1val, W2 = W2val, W3 = W3val, sA = seq(-4, 4, by = 0.2)))
  return(list(setWsA = datWset$sA, setWmat = setWmat))
}

get.density.sAdat <- function(nsamp = 100000) {
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
  # head(datO, 100)
}

def.nodeojb <- function(datO) {
  Kmax <- 1
  # nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"))
  nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"), nFnode = "nF")
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA")
  # netind_cl <- NetIndClass$new(Odata = datO)
  netind_cl <- NetIndClass$new(nobs = nrow(datO))
  # Define datNetObs:
  datnetW <- DatNet$new(netind_cl = netind_cl, nodes = nodes, VarNodes = nodes$Wnodes, addnFnode = TRUE)$make.sVar(Odata = datO, sVar.object = def_sW)
  datnetA <- DatNet$new(netind_cl = netind_cl, nodes = nodes, VarNodes = nodes$Anode)$make.sVar(Odata = datO, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
  return(list(datNetObs = datNetObs, netind_cl = netind_cl, def_sA = def_sA, def_sW = def_sW, nodes = nodes))
}

# directly fit a joint density for sA, sA2 ~ W, for sA - continuous
test.simple.fit.density.sA <- function() {
  nsamp <- 10000
  datO <- get.density.sAdat(nsamp)
  nodeobjs <- def.nodeojb(datO)
  testm.sW <- nodeobjs$def_sW$get.mat.sVar(data.df = datO, netind_cl = nodeobjs$netind_cl, addnFnode = "nF") # addnFnode = NULL)
  print("testm.sW"); print(head(testm.sW)); print("testm.sW map"); print(nodeobjs$def_sW$sVar.names.map)
  testm.sA <- nodeobjs$def_sA$get.mat.sVar(data.df = datO, netind_cl = nodeobjs$netind_cl)
  print("testm.sA"); print(head(testm.sA)); print("testm.sA map"); print(nodeobjs$def_sA$sVar.names.map)

  # Define est_params_list:
  reg.sVars <- list(outvars = c("sA"), predvars = c("W1", "W2", "W3"))
  # TO DO: Replace var with the index (column) such as, which(colnames(dat)%in%var); Then can use evaluator even when dat is matrix
  subsets_chr <- lapply(reg.sVars$outvars, function(var) {"!misfun("%+%var%+%")"})  # subsetting by !gvars$misval on sA:
  subsets_expr <- lapply(subsets_chr, function(subset_chr) { try(parse(text=subset_chr)[[1]]) }) # parses chr into a call
  sA_class <- nodeobjs$datNetObs$datnetA$type.sVar[reg.sVars$outvars]

  regclass <- RegressionClass$new(outvar.class = sA_class, outvar = reg.sVars$outvars, predvars = reg.sVars$predvars, subset = subsets_expr)
  summeas.g0 <- SummariesModel$new(reg = regclass, O.datnetA = nodeobjs$datNetObs$datnetA)
  summeas.g0$fit(data = nodeobjs$datNetObs)

  contsumobj <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  (intrvls <- contsumobj$intrvls)
  (intrvls.width <- diff(intrvls))
  length(intrvls.width)

  print(nodeobjs$datNetObs$get.sVar.intrvls("sA"))
  ord.sVar <- nodeobjs$datNetObs$discretize.sVar("sA")
  ord.sVar_bw <- intrvls.width[ord.sVar]
  print(head(cbind(sA = nodeobjs$datNetObs$dat.sVar[, "sA"], ord.sVar, bw = ord.sVar_bw, nodeobjs$datNetObs$dat.bin.sVar), 5))
  print("freq count for transformed ord.sVar: "); print(table(ord.sVar))
  plot(density(ord.sVar))
  hist(ord.sVar)

  summeas.g0$predict(newdata = nodeobjs$datNetObs)  # summeas.g0$sA_nms
  # Get P(sA|W) for the observed data (W,sA):
  # SHOULD BE SIMILAR TO THE OBSERVED DENSITY OF s.A (but discretized)
  h_gN <- summeas.g0$predictAeqa(obs.DatNet.sWsA = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  # plot(datO[,"sA"], h_gN, type = "p", cex = .5)
  # setWvals <- c(W1 = 0L, W2 = 0L, W3 = 0L)
  # subs000 <- (datO$W1==0L & datO$W2==0L & datO$W3==0L)
  setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
  # subs001 <- (datO$W1==0 & datO$W2==0 & datO$W3==1)
  # # setWvals <- c(W1 = 0, W2 = 1, W3 = 1)
  # subs011 <- (datO$W1==0 & datO$W2==1 & datO$W3==1)
  # setWvals <- c(W1 = 1, W2 = 1, W3 = 1)
  # subs111 <- (datO$W1==1 & datO$W2==1 & datO$W3==1)
  # setWvals <- c(W1 = 1, W2 = 0, W3 = 1)
  # subs100 <- (datO$W1==1 & datO$W2==0 & datO$W3==0)
  # setWvals <- c(W1 = 1, W2 = 0, W3 = 0)
  # subs100 <- (datO$W1==1 & datO$W2==0 & datO$W3==0)

  subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])
  sum(subs)
  setWdat_res <- get.setW.sAdat(setWvals, nsamp)

  length(unique(h_gN[subs]))
  sum(unique(h_gN[subs]))

  plot(datO[subs,"sA"], h_gN[subs], type = "p", cex = .3, col = "red")
  lines(density(setWdat_res$setWsA))

  plot(density(setWdat_res$setWsA))
  lines(datO[subs,"sA"], h_gN[subs], type = "p", cex = .3, col = "red")
  setWdat_res$setWmat

  # setWdat_res$setWmat
  # Get P(sA=sa|W) for fixed W and sa \in range(sA):
  # nodeobjs$datNetObs$mat.sVar <- as.matrix(setWdat_res$setWmat)
  # summeas.g0$predict(newdata = nodeobjs$datNetObs)
  # dens_p <- summeas.g0$predictAeqa(obs.DatNet.sWsA = nodeobjs$datNetObs)
  # plot(setWdat_res$setWmat[,"sA"], dens_p, type = "l")
}


# ---------------------------------------------------------------------------------------------------------
# Test iid TMLE fit for continous sA 
# ---------------------------------------------------------------------------------------------------------
# TMLE for causal effect in i.i.d. data with continuous exposure under continuous stochastic intervention;
# intervention g.star is defined by shifting the normal density of observed sA until g.star/g.0 >= 10, 
# then its truncated to be equal to g.0
get.iid.densityOdat <- function(nsamp = 100000, rndseed = NULL, trunc.const = 10, shift.const = 2) {
  library(simcausal)
  # nsamp = 100000
  # rndseed = 12345
  # trunc.const <- 10
  # shift.const <- 2
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("shift",  distr = "rconst", const = .(shift.const)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("r.obs.sA",  distr = "rconst", const = exp(shift * (sA - sA.mu - shift / 2))) +
      node("trunc.c",  distr = "rconst", const = .(trunc.const)) +
      node("untrunc.sA.gstar",  distr = "rconst", const = sA + shift) +
      # node("untrunc.sA.gstar",  distr = "rnorm", mean = sA.mu + shift, sd = 1) +
      # node("p.gstar.sA", distr = "rconst", const = (1/sqrt(2*.(pi))) * exp((-1/2) * (untrunc.sA.gstar - sA.mu)^2)) +
      # node("p.gstar.sA.gstar", distr = "rconst", const = (1/sqrt(2*.(pi))) * exp((-1/2) * (untrunc.sA.gstar - (sA.mu + shift))^2)) +
      node("r.new.sA",  distr = "rconst", const = exp(shift * (untrunc.sA.gstar - sA.mu - shift / 2))) +
      node("trunc.sA.gstar",  distr = "rconst", const = ifelse(r.new.sA > trunc.c, sA, untrunc.sA.gstar)) +
      node("probY", distr = "rconst", const = plogis(-0.45 * sA - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y", distr = "rbern", prob = plogis(-0.45 * sA - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) + 
      node("probY.gstar", distr = "rconst", const = plogis(-0.45 * trunc.sA.gstar - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y.gstar", distr = "rbern", prob = plogis(-0.45 * trunc.sA.gstar - 0.5 * W1 - 0.58 * W2 - 0.33 * W3))

  D <- set.DAG(D)
  datO <- sim(D, n = nsamp, rndseed = rndseed)

  head(datO, 100)
  print("mean(datO$Y): " %+% mean(datO$Y));
  # [1] 0.31625
  psi0 <- mean(datO$Y.gstar)
  print("psi0: " %+% psi0)
  # [1] 0.22378

  return(list(psi0 = psi0, datO = datO))
}

run.1sim.tmlenet <- function(nsamp, Qform, f.gstar, trunc.const = 10, shift.const = 2) {
  datO <- get.iid.densityOdat(nsamp = nsamp, rndseed = NULL, trunc.const = trunc.const, shift.const = shift.const)$datO
  datO <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  Kmax <- 1
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA")

  tmlenet_res <- tmlenet(data = datO, Anode = "sA", Wnodes = c("W1", "W2", "W3"), Ynode = "Y",
                          Kmax = Kmax, 
                          # IDnode = "ID", NETIDnode = "ID",
                          nFnode = NULL,
                          f_gstar1 = f.gstar,
                          sW = def_sW, sA = def_sA,
                          Qform = Qform,
                          # correct Q:
                          # Qform = "Y ~ W1 + W2 + W3 + sA",
                          # misspecified Q:
                          # Qform = "Y ~ W3 + sA",

                          hform = "sA ~ W1 + W2 + W3",
                          hform.gstar = "sA ~ W1 + W2 + W3",
                          gform = "sA ~ W1 + W2 + W3",

                          optPars = list(
                            # f_g0 = f.A, # tested, works
                            n_MCsims = 10
                          )
                          )

  CIs <- tmlenet_res$EY_gstar1$CIs
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

# Run one TMLE simulation for iid data sampled from get.iid.densityOdat, estimating psi0 under trunced g.star
# Used for Verifying: 1) consistency, 2) double robustness, 3) correct as. coverage
test.onesim.iid.tmlefit <- function() {
  trunc.const <- 10
  shift.const <- 2

  # get true psi.0:
  datO <- get.iid.densityOdat(nsamp = 2000000, rndseed = 12345, trunc.const = trunc.const, shift.const = shift.const)
  psi0 <- datO$psi0
  print("psi0: " %+% psi0)
  # [1] "psi0: 0.2241085"

  nsamp <- 2000
  # correct Q:
  Qform.corr <- "Y ~ W1 + W2 + W3 + sA"
  # misspecified Q:
  Qform.mis <- "Y ~ W3 + sA"
  Qforms <- c(Qform.corr = Qform.corr, Qform.mis = Qform.mis)

  f.gstar <- function(data, ...) {
    sA.mu <- 0.98 * data[,"W1"] + 0.58 * data[,"W2"] + 0.33 * data[,"W3"]
    untrunc.sA <- rnorm(n = nrow(data), mean = sA.mu + shift.const, sd = 1)
    r.new.sA <- exp(shift.const * (untrunc.sA - sA.mu - shift.const / 2))
    trunc.sA <- ifelse(r.new.sA > trunc.const, untrunc.sA - shift.const, untrunc.sA)
    return(trunc.sA)
  }

  res <- run.1sim.tmlenet(nsamp = nsamp, Qform = Qform.mis, f.gstar = f.gstar,  trunc.const = 10, shift.const = 2)
  res
}