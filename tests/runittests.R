### --- Test setup ---

if(FALSE) {
  # DEBUGGING:
  # memory & speed profiling
  # Rprof("tmle_run_memuse.out", memory.profiling=TRUE)
  library("RUnit")
  library("roxygen2")
  library("devtools")
  setwd(".."); getwd()
  document()
  setwd("..")
  install("tmlenet", build_vignettes=FALSE) # INSTALL W/ devtools:

  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf tmlenet")  # just create the pdf manual from help files

  # CHECK AND BUILD PACKAGE:
  # setwd(".."); setwd(".."); getwd()
  # devtools::check() # runs check with devtools

  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  devtools::build(args = "--compact-vignettes") # build package tarball compacting vignettes
  # devtools::build() # build package tarball
  setwd("..")
  system("R CMD check --as-cran tmlenet_0.0.9.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes tmlenet") # check without building the pdf manual and not building vignettes
      ## system("R CMD build tmlenet --no-build-vignettes")
      ## system("R CMD build tmlenet")
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./tmlenet_0.0.9.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(tmlenet)
  # tmlenet:::addvectorfcn("poisson")
  # tmlenet:::debug_set() # SET TO DEBUG MODE
  # tmlenet:::debug_off() # SET DEBUG MODE OFF

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

psi_RDs_DAG2a <- NULL
psi_RDs_DAG2b <- NULL

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

`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))


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
  netind_cl <- NetIndClass$new(Odata = datO)
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

  summeas.g0 <- SummariesModel$new(sA_class = sA_class, # (1) TO BE REPLACED WITH RegressionClass object that defines sA_class, sA_nms, sW_nms
                                    sA_nms = reg.sVars$outvars,
                                    sW_nms = reg.sVars$predvars,
                                    subset = subsets_expr,
                                    O.datnetA = nodeobjs$datNetObs$datnetA)

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

  # Algorithm for obtaining the normalization constants for unequal bin lengths (bin_bymass = TRUE)
  # NOT NEEDED: we are estimating the discretized density, not the pmf, hence no need to normalize when P(bin_i) * (1/bini_width)
    # 0) New class BinDatWt that inherits from BinDat that calculates the weights below
    # 1) for each bin B_i calculate p1_B.i_full = P(B_i = 1 | X_mat) for ALL observations (nrow(X_mat) = nobs)
    # 2) for each bin B_i calculate the joint p1_ord.sVar_B.i = P(sVar in B_i | X_mat) = prod_{j=1,...,i-1}(1 - p1_B.j_full) * p1_B.i_full
    # 3) for each h_i (bin B_i length), take a sum (1 / h_i) * p1_ord.sVar_B.i, obtaining an n-vector of normalization constants
  # allPsAsW.models <- contsumobj$getPsAsW.models()
  # # allPsAsW.models$`P(sA|sW).1`$getfit
  # # sum(allPsAsW.models$`P(sA|sW).1`$bindat$subset_idx)
  # # length(allPsAsW.models$`P(sA|sW).1`$bindat$getY)
  # # head(allPsAsW.models$`P(sA|sW).1`$getprobA1)
  # p.ord.sVar_B.i_list <- list()
  # p0.prev_B.i <- matrix(1L, nrow = nrow(datO), ncol = 1)
  # for (PsAsWmodel in allPsAsW.models) {  
  #   coefs <- PsAsWmodel$getfit$coef
  #   predvars <- PsAsWmodel$reg$predvars
  #   newdata <- nodeobjs$datNetObs
  #   Xmat_all <- cbind(1, newdata$get.dat.sWsA(covars = predvars))
  #   p1.B.i <- logitlinkinv(Xmat_all[,!is.na(coefs), drop = FALSE] %*% coefs[!is.na(coefs)])
  #   if (all(is.na(coefs))) p1.B.i[] <- 1L
  #   p.ord.sVar_B.i <- p1.B.i * p0.prev_B.i
  #   p0.prev_B.i <- p0.prev_B.i * (1L - p1.B.i)
  #   p.ord.sVar_B.i_list <- c(p.ord.sVar_B.i_list, list(p.ord.sVar_B.i))
  # }
  # p.ord.sVar_mat <- do.call("cbind", p.ord.sVar_B.i_list)
  # p.ord.sVar_bw_mat <- p.ord.sVar_mat / matrix(intrvls.width, nrow = nrow(p.ord.sVar_mat), ncol = ncol(p.ord.sVar_mat), byrow = TRUE)
  # p.ord.sVar.normconsts <- rowSums(p.ord.sVar_bw_mat)
}

get.iid.densityOdat <- function(nsamp = 100000, rndseed = NULL, trunc.const = 10, shift.const = 2) {
  `%+%` <- function(a, b) paste0(a, b)
  library(simcausal)
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

  # summary(datO$untrunc.sA.gstar)
  # summary(datO$trunc.sA.gstar)
  # sum(datO$r.new.sA > trunc.const)
  # # sA under g0 vs. g.star:
  # plot(density(datO$sA), xlim = c(-5, 8), ylim = c(0,1)) # observed g.0
  # lines(density(datO$untrunc.sA.gstar), col = "red", lty="solid") # untruncated g.star
  # lines(density(datO$trunc.sA.gstar), col = "red", lty="dashed") # truncated g.star

  # # Distribution of observed probY by sA cut-offs:
  # summary(datO$probY)
  # plot(density(datO$probY), ylim = c(0, 8))
  # lines(density(datO$probY[datO$sA>-1]), col = "red", lty="solid")
  # lines(density(datO$probY[datO$sA>0]), col = "red", lty="dashed")
  # lines(density(datO$probY[datO$sA>1]), col = "red", lty="dotted")
  # lines(density(datO$probY[datO$sA>2]), col = "red", lty="dotdash")

  # # Distribution of probY under g.star by sA cut-offs:
  # summary(datO$probY.gstar)
  # plot(density(datO$probY.gstar), ylim = c(0, 8))
  # lines(density(datO$probY.gstar[datO$sA>-1]), col = "red", lty="solid")
  # lines(density(datO$probY.gstar[datO$sA>0]), col = "red", lty="dashed")
  # lines(density(datO$probY.gstar[datO$sA>1]), col = "red", lty="dotted")
  # lines(density(datO$probY.gstar[datO$sA>2]), col = "red", lty="dotdash")

  # # Effect of W1 on sA:
  # plot(density(datO$sA))
  # lines(density(datO$sA[datO$W1==0]), col = "red")
  # lines(density(datO$sA[datO$W1==1]), col = "red", lty="dashed")
  # # Effect of W1 on probY:
  # plot(density(datO$probY))
  # lines(density(datO$probY[datO$W1==0]), col = "red")
  # lines(density(datO$probY[datO$W1==1]), col = "red", lty="dashed")
  # # Effect of W2 on sA:
  # plot(density(datO$sA))
  # lines(density(datO$sA[datO$W2==0]), col = "red")
  # lines(density(datO$sA[datO$W2==1]), col = "red", lty="dashed")
  # # Effect of W2 on probY:
  # plot(density(datO$probY))
  # lines(density(datO$probY[datO$W2==0]), col = "red")
  # lines(density(datO$probY[datO$W2==1]), col = "red", lty="dashed")
  # # Effect of W3 on sA:
  # plot(density(datO$sA))
  # lines(density(datO$sA[datO$W3==0]), col = "red")
  # lines(density(datO$sA[datO$W3==1]), col = "red", lty="dashed")
  # # Effect of W3 on probY:
  # plot(density(datO$probY))
  # lines(density(datO$probY[datO$W3==0]), col = "red")
  # lines(density(datO$probY[datO$W3==1]), col = "red", lty="dashed")

  return(list(psi0 = psi0, datO = datO))
}


run.sim.tmlenet <- function(nsamp, Qform, f.gstar) {
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
                          Qform.new = Qform,
                          # correct Q:
                          # Qform.new = "Y ~ W1 + W2 + W3 + sA",
                          # misspecified Q:
                          # Qform.new = "Y ~ W3 + sA",

                          hform.new = "sA ~ W1 + W2 + W3",
                          hform.gstar.new = "sA ~ W1 + W2 + W3",
                          gform.new = "sA ~ W1 + W2 + W3",

                          opt.params = list(
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

# Run TMLE for iid data sampled from get.iid.densityOdat, estimating psi0 under trunced g.star
# Used for Verifying: 1) consistency, 2) double robustness, 3) correct as. coverage
test.iid.tmlefit <- function() {
  `%+%` <- function(a, b) paste0(a, b)
  library(devtools)
  load_all("../") # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  trunc.const <- 10
  shift.const <- 2

  # get true psi.0:
  datO <- get.iid.densityOdat(nsamp = 2000000, rndseed = 12345, trunc.const = trunc.const, shift.const = shift.const)
  psi0 <- datO$psi0
  print("psi0: " %+% psi0)
  # [1] "psi0: 0.2241085"

  nsamp <- 2000
  # correct Q:
  Qform.corr = "Y ~ W1 + W2 + W3 + sA"
  # misspecified Q:
  Qform.mis = "Y ~ W3 + sA"
  Qforms <- c(Qform.corr = Qform.corr, Qform.mis = Qform.mis)

  f.gstar <- function(data, ...) {
    sA.mu <- 0.98 * data[,"W1"] + 0.58 * data[,"W2"] + 0.33 * data[,"W3"]
    untrunc.sA <- rnorm(n = nrow(data), mean = sA.mu + shift.const, sd = 1)
    r.new.sA <- exp(shift.const * (untrunc.sA - sA.mu - shift.const / 2))
    trunc.sA <- ifelse(r.new.sA > trunc.const, untrunc.sA - shift.const, untrunc.sA)
    return(trunc.sA)
  }
  res <- run.sim.tmlenet(nsamp = nsamp, Qform = Qform.mis, f.gstar = f.gstar)


  source("determineParallelBackend.R")
  require(doParallel)
  registerDoParallel(cores=1)
  nsim <- 50 # nsim <- 500

  runsim <- function(nsim, scen.name) {
    Qform <- Qforms[scen.name]
    print("Qform: " %+% Qform)
    simRes <- foreach(t.counter = 1 : nsim) %dopar% {
      run.sim.tmlenet(nsamp = nsamp, Qform = Qform, f.gstar = f.gstar)
    }
    save(list=c("simRes"), file="sims." %+% scen.name %+% ".RData")
    return(simRes)
  }
  simRes_all <- list()
  for (scen.name in names(Qforms)) {
    # scen.name <- names(Qforms)[2]
    simres <- runsim(nsim = nsim, scen.name = scen.name)
    (ests_sims <- sapply(simres, '[[', "est"))
    abs.bias <- abs(rowMeans(ests_sims) - psi0)[c("tmle_B", "h_iptw", "mle")]
    cover_sims <- sapply(simres, '[[', "cover")
    coverage <- rowMeans(cover_sims)[c("tmle_B", "h_iptw")]
    simRes <- list(scen.name = scen.name, psi0 = psi0, nsim = nsim, nsamp = nsamp, abs.bias = abs.bias, coverage = coverage)
    simRes_all <- c(simRes_all, simRes)
    print("simRes_all"); print(simRes_all)
  }
  save(list=c("simRes_all"), file="res.sims.nsims"%+%nsim%+%".RData")

}




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
                          Qform.new = "Y2_1 ~ W1 + W2 + W3 + A1_0 + A2_0 + A1_1 + A2_1 + Y1_1",
                          hform.new = "A1_1 + A2_1 + Y1_1 ~ W1 + W2 + W3 + A1_0 + A2_0",
                          hform.gstar.new = "A1_1 + A2_1 + Y1_1 ~ W1 + W2 + W3 + A1_0 + A2_0",
                          gform.new = "A2 ~ W1 + W2 + W3 + A1_0 + A2_0",
                          opt.params = list(
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
# TESTING ContinSummaryModel class and NewSummaryModel.contin constructor
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


test.RegressionClass <- function() {
# TEST RegressionClass:
reg_test <- RegressionClass$new(outvar.class = gvars$sVartypes$bin, 
                                  outvar = "A", 
                                  predvars = c("W1", "W2"), 
                                  subset = expression(A=0))
  # class(reg_test)
  reg_test$reg
}


test.continous.sA <- function() {
    # I) Build network vectors: (W, W_netF_1, ..., W_netF_k) for each W in Wnodes by PRE-ALLOCATING netW_full:
    datnetW <- DatNet$new(Odata = data, NetInd_k = NetInd_k, Kmax = k, nodes = node_l, VarNodes = node_l$Wnodes, addnFnode = TRUE)
    # II) APPLY THE SUMMARY MEASURE FUNCTIONS / EXPRESSION TO netW_full to OBTAIN sW columns SELECT ONLY sW columns in hform_g0 and hfrom_gstar or use all? 
    obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms, norm.c.sVars = TRUE)$dat.sVar
    datnetW$def_cbin_intrvls()
    print("Detected intervals: "); print(datnetW$cbin_intrvls)
    print("Detected nbins: "); print(datnetW$nbins)
    oldncats1 <- set.maxncats(5)
    oldnbins1 <- set.nbins(10)
    print("No normalization. Binning by mass")
      obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms)$dat.sVar
      print("head(obsdat.sW)"); print(head(obsdat.sW))
      print("Var types: "); str(datnetW$type.sVar)
      print("Contin covars names: "); print(datnetW$names.c.sVar)
      defints1 <- datnetW$def_cbin_intrvls()$cbin_intrvls
      print("No normalization bin intervals by mass: "); print(defints1)
      print("nbins: "); print(datnetW$nbins); 
      print("Testing ordinals with ncats < nbins get nbins = ncats:"); print(datnetW$nbins["nFriends"] < gvars$nbins)

    print("No normalization. Binning by equal length")
      intlen_ints <- datnetW$def_cbin_intrvls(bin_bymass = FALSE)$cbin_intrvls
      print("No normalization bin intervals by length: "); print(intlen_ints)
      print("nbins: "); print(datnetW$nbins)

    print("Testing ordinals with ncats > nbins get collapsed into fewer cats:")
    set.nbins(4)
      obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms)$dat.sVar
      defints2 <- datnetW$def_cbin_intrvls()$cbin_intrvls
      print("New bins with collapsed ordinals: "); print(defints2)
      print("nbins: "); print(datnetW$nbins)

    set.nbins(10)
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
      newint <- seq(from = 0, to = 1, length.out = (gvars$nbins+1))
      datnetW$def_cbin_intrvls(cbin_intrvls = newint)
      print("Assigned bin intervals 1: "); print(datnetW$cbin_intrvls)
      print("nbins: "); print(datnetW$nbins)
      print("Correctly assigned bin intervals 1: "); print(all.equal(as.vector(datnetW$cbin_intrvls[[1]]), as.vector(datnetW$cbin_intrvls[[2]])))

    print("testing overwriting of bin intervals 2:")
      int1 <-  seq(from = 0, to = 0.5, length.out = (gvars$nbins+1))
      int2 <-  seq(from = 0.5, to = 1, length.out = (gvars$nbins+1))
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
    set.maxncats(oldncats1)
    set.nbins(oldnbins1)
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


test.faster_tolongdata <- function() {
    #-------------------------------------------------------------
    # Testing the converter to long format that is based on package data.table
    #-------------------------------------------------------------
    # t_end <- 16
    # library(simcausal)

    # D <- DAG.empty()
    # D <- D + node("L2", t=0,        distr="rbern", prob=0.05, order=1)
    # D <- D + node("L1", t=0,        distr="rbern", prob=ifelse(L2[0]==1,0.5,0.1), order=2)
    # D <- D + node("A1", t=0,        distr="rbern", prob=ifelse(L1[0]==1 & L2[0]==0, 0.5, ifelse(L1[0]==0 & L2[0]==0, 0.1, ifelse(L1[0]==1 & L2[0]==1, 0.9, 0.5))), order=3)
    # D <- D + node("A2", t=0,        distr="rbern", prob=0, order=4, EFU=TRUE)
    # D <- D + node("Y",  t=0,        distr="rbern", prob=plogis(-6.5 + L1[0] + 4*L2[0] + 0.05*I(L2[0]==0)), order=5, EFU=TRUE)
    # D <- D + node("L2", t=1:t_end,  distr="rbern", prob=ifelse(A1[t-1]==1, 0.1, ifelse(L2[t-1]==1, 0.9, min(1,0.1 + t/16))), order=6+4*(0:(t_end-1)))
    # D <- D + node("A1", t=1:t_end,  distr="rbern", prob=ifelse(A1[t-1]==1, 1, ifelse(L1[0]==1 & L2[0]==0, 0.3, ifelse(L1[0]==0 & L2[0]==0, 0.1, ifelse(L1[0]==1 & L2[0]==1, 0.7, 0.5)))), order=7+4*(0:(t_end-1)))
    # D <- D + node("A2", t=1:t_end,  distr="rbern", prob=plogis(-3.5 + 0.5*A1[t]+0.5*L2[t]), order=8+4*(0:(t_end-1)), EFU=TRUE) # informative censoring

    # # this takes longer (6 sec longer for 1Mil obs)
    # # D <- D + node("Y",  t=1:t_end,  distr="rbern", prob=plogis(-6.5 + L1[0] + 4*L2[t] + 0.05*sum(I(L2[0:t]==rep(0,(t+1))))), order=9+4*(0:(t_end-1)), EFU=TRUE)
    # D <- D + node("Y",  t=1:t_end,  distr="rbern", prob=plogis(-6.5 + L1[0] + 4*L2[t] + 0.05*(sum(L2[0:t])==0)), order=9+4*(0:(t_end-1)), EFU=TRUE)
    # lDAG3 <- set.DAG(D)

    # #-------------------------------------------------------------
    # # Adding dynamic actions (indexed by a real-valued parameter)
    # #-------------------------------------------------------------
    # # act_t0_theta <- node("A1",t=0, distr="rbern", prob=ifelse(L2[0] >= theta,1,0))
    # act_t0_theta <- node("A1",t=0, distr="rbern", prob=ifelse(L2[0] >= theta,1,0))
    # act_tp_theta <- node("A1",t=1:t_end, distr="rbern", prob=ifelse(A1[t-1]==1,1,ifelse(L2[t] >= theta,1,0)))
    # act_NoCens <- node("A2",t=0:t_end, distr="rbern", prob=0)
    # actionnodes <- c(act_t0_theta, act_tp_theta, act_NoCens)
    # D <- lDAG3 + action("A1_th0", nodes=actionnodes, theta=0)
    # D <- D + action("A1_th1", nodes=actionnodes, theta=1)

    # #-------------------------------------------------------------
    # # Testing conversion of observed data to long format 
    # #-------------------------------------------------------------
    # # NO carry forward imputation:
    # O_dat_df <- simobs(D, n=500, rndseed = 123)
    # system.time(O_dat_long <- DF.to.long(O_dat_df))
    # system.time(O_dat_long_DT <- DF.to.longDT(O_dat_df))

    # checkIdentical(O_dat_long$ID, O_dat_long_DT$ID)
    # checkIdentical(O_dat_long$L1_0, O_dat_long_DT$L1_0)
    # checkIdentical(O_dat_long$t, O_dat_long_DT$t)
    # checkIdentical(O_dat_long$L2, O_dat_long_DT$L2)
    # checkIdentical(O_dat_long$A1, O_dat_long_DT$A1)
    # checkIdentical(O_dat_long$A2, O_dat_long_DT$A2)
    # checkIdentical(O_dat_long$Y, O_dat_long_DT$Y)

    # # With carry forward imputation of Y (vs 1):
    # O_dat_df <- simobs(D, n=500, rndseed = 123)
    # O_dat_LTCF <- doLTCF(data=O_dat_df, LTCF="Y")
    # system.time(O_dat_long_LTCF_v1 <- DF.to.long(O_dat_LTCF))
    # system.time(O_dat_long_DT_LTCF_v1 <- DF.to.longDT(O_dat_LTCF))

    # checkIdentical(O_dat_long_LTCF_v1$ID, O_dat_long_DT_LTCF_v1$ID)
    # checkIdentical(O_dat_long_LTCF_v1$L1_0, O_dat_long_DT_LTCF_v1$L1_0)
    # checkIdentical(O_dat_long_LTCF_v1$t, O_dat_long_DT_LTCF_v1$t)
    # checkIdentical(O_dat_long_LTCF_v1$L2, O_dat_long_DT_LTCF_v1$L2)
    # checkIdentical(O_dat_long_LTCF_v1$A1, O_dat_long_DT_LTCF_v1$A1)
    # checkIdentical(O_dat_long_LTCF_v1$A2, O_dat_long_DT_LTCF_v1$A2)
    # checkIdentical(O_dat_long_LTCF_v1$Y, O_dat_long_DT_LTCF_v1$Y)

    # # With carry forward imputation of Y (vs 2):
    # O_dat_df_LTCF <- simobs(D, n=500, LTCF="Y", rndseed = 123)
    # system.time(O_dat_long_LTCF <- DF.to.long(O_dat_df_LTCF))
    # system.time(O_dat_long_DT_LTCF <- DF.to.longDT(O_dat_df_LTCF))

    # checkIdentical(O_dat_long_LTCF$ID, O_dat_long_DT_LTCF$ID)
    # checkIdentical(O_dat_long_LTCF$L1_0, O_dat_long_DT_LTCF$L1_0)
    # checkIdentical(O_dat_long_LTCF$t, O_dat_long_DT_LTCF$t)
    # checkIdentical(O_dat_long_LTCF$L2, O_dat_long_DT_LTCF$L2)
    # checkIdentical(O_dat_long_LTCF$A1, O_dat_long_DT_LTCF$A1)
    # checkIdentical(O_dat_long_LTCF$A2, O_dat_long_DT_LTCF$A2)
    # checkIdentical(O_dat_long_LTCF$Y, O_dat_long_DT_LTCF$Y)

    # #-------------------------------------------------------------
    # # Testing conversion of full data to long format (with carry forward imputation)
    # #-------------------------------------------------------------
    # X_dat <- simfull(A(D), n=500, rndseed = 123)
    # attributes(X_dat[[1]])$node_nms

    # system.time(X_dat_l <- lapply(X_dat, DF.to.long))
    # system.time(X_dat_lDT <- lapply(X_dat, DF.to.longDT))

    # checkIdentical(X_dat_l[["A1_th0"]]$ID, X_dat_lDT[["A1_th0"]]$ID)
    # checkIdentical(X_dat_l[["A1_th0"]]$L1_0, X_dat_lDT[["A1_th0"]]$L1_0)
    # checkIdentical(X_dat_l[["A1_th0"]]$t, X_dat_lDT[["A1_th0"]]$t)
    # checkIdentical(X_dat_l[["A1_th0"]]$L2, X_dat_lDT[["A1_th0"]]$L2)
    # checkIdentical(X_dat_l[["A1_th0"]]$A1, X_dat_lDT[["A1_th0"]]$A1)
    # checkIdentical(X_dat_l[["A1_th0"]]$A2, X_dat_lDT[["A1_th0"]]$A2)
    # checkIdentical(X_dat_l[["A1_th0"]]$Y, X_dat_lDT[["A1_th0"]]$Y)

    # checkIdentical(X_dat_l[["A1_th1"]]$ID, X_dat_lDT[["A1_th1"]]$ID)
    # checkIdentical(X_dat_l[["A1_th1"]]$L1_0, X_dat_lDT[["A1_th1"]]$L1_0)
    # checkIdentical(X_dat_l[["A1_th1"]]$t, X_dat_lDT[["A1_th1"]]$t)
    # checkIdentical(X_dat_l[["A1_th1"]]$L2, X_dat_lDT[["A1_th1"]]$L2)
    # checkIdentical(X_dat_l[["A1_th1"]]$A1, X_dat_lDT[["A1_th1"]]$A1)
    # checkIdentical(X_dat_l[["A1_th1"]]$A2, X_dat_lDT[["A1_th1"]]$A2)
    # checkIdentical(X_dat_l[["A1_th1"]]$Y, X_dat_lDT[["A1_th1"]]$Y)

    # BENCHMARKING for 50K: gain of x5.4 factor
    # old convert to long
    # user  system elapsed 
    # 13.943   1.783  15.766 
    # new DF.to.longDT
    # user  system elapsed 
    # 2.564   0.370   2.935 

    # BENCHMARKING for 500K: gain of x5.4 factor
    # old convert to long
    # user  system elapsed 
    # 140.378  18.398 158.853     
    # new DF.to.longDT    
    # user  system elapsed 
    # 28.753   4.092  32.844 

    # CONVERTING BACK TO WIDE FORMAT:
    # ## convert long to wide using dcast from data.table
    # dcast.data.table(data, formula, fun.aggregate = NULL, 
    #     ..., margins = NULL, subset = NULL, fill = NULL, 
    #     drop = TRUE, value.var = guess(data),
    #     verbose = getOption("datatable.verbose"))


    #-------------------------------------------------------------
    # BENCHMARKING current package rowSums with data.table version - abandoned for now
    #-------------------------------------------------------------
    # t_pts <- 0:16
    # t_vec <- "_"%+%(t_pts)
    # (L2_names <- "L2"%+%t_vec)
    # # library(data.table)
    # system.time( X_dat_1 <- simfull(A(D)[1], n=1000000, LTCF="Y", rndseed = 123)[[1]])
    # X_dat_1 <- X_dat_1[,-1]
    # nrow(X_dat_1)
    # colnames(X_dat_1)
    # head(X_dat_1)
    # X_dat_1_DT = data.table(X_dat_1)
    # setkey(X_dat_1_DT,ID)
    # L2_idx <- which(names(X_dat_1_DT)%in%L2_names)
    # ID_idx <- which(names(X_dat_1_DT)%in%"ID")
    # ncol(X_dat_1_DT)
    # # COMPARING data.table to current data.rame rowSums
    # # fast version 1  of row sums 1
    # system.time(
    # X_dat_1_DT[, SumRow := rowSums(.SD), .SDcols = L2_idx]
    # )
    # user  system elapsed 
    # 0.139   0.030   0.181     

    # # faster version 2 of row sums
    # system.time(
    # X_dat_1_DT[, SumRow2 := Reduce(`+`, .SD), .SDcol = L2_idx]
    # )
    #  user  system elapsed 
    # 0.075   0.027   0.120

    # # version 3 using set
    # # i  In set(), integer row numbers to be assigned value.
    # # NULL represents all rows more efficiently than creating a vector such as 1:nrow(x).
    # # j In set(), integer column number to be assigned value.    
    # # x - DT, i - row indx, j - col indx, value - new val
    # set(x, i=NULL, j, value)
    # system.time(for (i in 1:nrow(X_dat_1_DT)) set(X_dat_1_DT,as.integer(i),"SumRow3",i))

    # # current version of row sums (slowest)
    # system.time(newVar_vold <- rowSums(X_dat_1[,L2_names]))
    # # COMPARING data.table to current data.rame rowSums with some column operations
    # system.time(X_dat_1_DT[, SumRow := rowSums(.SD==c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), .SDcols = L2_idx])
    # system.time(newVar_vold <- rowSums(I(X_dat_1[,L2_names] == cbind(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))))

}

