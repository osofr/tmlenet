#----------------------------------------------------------------------------------
# Estimate IC-based Variance (as. Var) and CIs (based on f_W and fY)
# New fast method for as var calculation (matrix vs)
#----------------------------------------------------------------------------------
# use estimates of fWi (hold all W's fixed at once),
# loop over all intersecting friends networks
# calculate R_ij*(fW_i*fW_j) - see page 33 Network paper vdL
#----------------------------------------------------------------------------------
# unified estimator naming used throughout the package:
# c("TMLE", "h_IPTW", "MLE")
#----------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# THIS FUNCTION IS NOT USED, ONLY HERE FOR TESTING SPARSE CONNECTIVITY MATRIX CONVERSION / CREATION
get_sparse_Fiintersectmtx <- function() {
  # ------------------------------------------------------------------------------
  # Simulate some network data
  # ------------------------------------------------------------------------------
  n <- 20000
  # n <- 10000
  # n <- 1000
  # fvec_i <- abs(rnorm(n))
  fvec_i <- rnorm(n)
  # fvec_i <- rep_len(1,n)
  # NetInd_k <- matrix(NA, nrow=n, ncol=3)
  # NetInd_k[1,] <- c(2,3,NA)
  # NetInd_k[4,] <- c(2,NA,NA)
  # nF <- rep.int(0,n)
  # nF[1] <- 2; nF[4] <- 1
  # Kmax <- 150
  Kmax <- 200
  # Kmax <- 15
  NetInd_k <- t(replicate(n, sample(1:n, Kmax, replace = FALSE)))
  nF <- rep.int(Kmax, n)
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Using sparse matrix implementation to get conn_ind_mtx_1st_indir
  # ------------------------------------------------------------------------------
  # estvartime <- system.time({
    # NetInd as sparse adjacency matrix (new version returns pattern sparse mat nngCMatrix):
    sparse_mat <- NetInd.to.sparseAdjMat(NetInd_k, nF = nF, add_diag = TRUE)
    # Second pass over columns of connectivity mtx to connect indirect intersections (i and j have a common friend but are not friends themselves):
    sparse_mat_ind <- Matrix::crossprod(sparse_mat) # t(sparse_mat)%*%sparse_mat returns nsCMatrix (only non-zero entries)
    # MOVED TO est.sigma_sparse():
    Dstar_crossprod_new <- sum_crossprod_Fij(sparseAdjMat = sparse_mat_ind, fvec_i = fvec_i)
    Dstar_new <- (1/n) * ((2*Dstar_crossprod_new) + sum(fvec_i^2))
  # })
  # print("estvartime"); print(estvartime)

  # print(object.size(NetInd_k), units="Mb")
  # 3906.4 Kb # 0 Gb
  # sparse_mat as dgCMatrix:
  # print(object.size(sparse_mat), units="Mb")
  # 12,032.1 Kb
  # system.time({
  #   conn_ind_mtx_1st_indir <- get.Fiintersectmtx(n = n)$conn_ind_mtx_1st
  #   Dstar_old <- est.sigma_fsum(get.crossprodmtx(fvec_i), conn_ind_mtx_1st_indir)        
  # })
  # print(all.equal(Dstar_new, Dstar_old))
  # # print(all.equal(Dstar, Dstar_old)); print(Dstar);
  # print(Dstar_new); print(Dstar_old)
  return(list(conn_ind_mtx_1st = sparse_mat_ind))
}

# sparse matrix class from Matrix package:
# dgCMatrix
  # class is a class of sparse numeric matrices in the compressed, sparse, column-oriented format.
# nsparseMatrix-classes
# ?'nsCMatrix-class'
# ngCMatrix, nsCMatrix, and ntCMatrix
  # class of sparse pattern matrices, i.e., binary matrices conceptually with TRUE/FALSE entries. Only the positions of the elements that are TRUE are stored.
  # Objects can be created by calls of the form new("ngCMatrix", ...) and so on. 
  # More frequently objects are created by coercion of a numeric sparse matrix to the pattern form for use in the symbolic 
  # analysis phase of an algorithm involving sparse matrices. Such algorithms often involve two phases: a symbolic phase wherein the 
  # positions of the non-zeros in the result are determined and a numeric phase wherein the actual results are calculated. 
  # During the symbolic phase only the positions of the non-zero elements in any operands are of interest, hence numeric sparse matrices can be treated as sparse pattern matrices.
# lsparseMatrix-classes {Matrix}
# lsCMatrix
# The second letter in the name of these non-virtual classes indicates general, symmetric, or triangular.
# The lsparseMatrix class is a virtual class of sparse matrices with TRUE/FALSE or NA entries. Only the positions of the elements that are TRUE are stored.

# Copied from simcausal, then modified to work as sparse pattern matrix:
# return pattern sparse matrix, no @x is recorded (ngMatrix):
NetInd.to.sparseAdjMat <- function(NetInd_k, nF, add_diag = FALSE) {
  nobs <- nrow(NetInd_k)
  sdims <- c(nobs, nobs)
  nnonzero <- sum(!is.na(NetInd_k))
  sparse_p <- as.integer(c(0, cumsum(nF)))
  sparse_iwNA <- as.vector(t(NetInd_k))
  sparse_i <- sparse_iwNA[!is.na(sparse_iwNA)] - 1
  out_sparseAdjMat <-  Matrix::sparseMatrix(i = sparse_i, p = sparse_p, dims = sdims, index1 = FALSE, giveCsparse = TRUE)
  if (add_diag) diag(out_sparseAdjMat) <- TRUE
  return(out_sparseAdjMat)
}

# For correlated i,j (s.t. i<j), add a term f_vec[i]*f_vec[j] to the cumulative sum
# Helper function to calculate cross product sum of correlated f_Wi (see p.33 of vdL)
sum_crossprod_Fij <- function(sparseAdjMat, fvec_i) {
  # sparseAdjMat:
    # nnzero: number of non-zero elements
    # i: These are the 0-based row numbers for each non-zero element in the matrix.
    # "integer" of length nnzero, 0- based row numbers for each non-zero element in the matrix, i.e., i must be in 0:(nrow(.)-1).
    # p: integer vector for providing pointers, one for each column, to the initial (zero-based) index of elements in the column.
    # .@p is of length ncol(.) + 1, with p[1] == 0 and
    # p[length(p)] == nnzero, such that in fact, diff(.@p) are the number of non-zero elements for each column.  
  assertthat::assert_that(is(sparseAdjMat, "sparseMatrix"))
  # 1) The number of friends for each observation:
  nF <- as.integer(diff(sparseAdjMat@p))
  # 2) Column based cummulative number of non-zero entries (cummulative nF)
  cumFindx <- sparseAdjMat@p
  # 3) All non-zero elements as a vector of 0-based row numbers:
  base0_IDrownums <- sparseAdjMat@i
  # 4) For each observation i that has non-zero nF (friends), add fvec_i[i]*fvec_Fj for each friend Fj of i:  
  non0nF.idx <- which(nF > 1L) # don`t care if nF[i]=1 since it means i has 0 actual friends (i itself is included in nF)
  # non0nF.idx <- which(nF > 0L)
  Dstar_crossprod <- 0
  for (idx in non0nF.idx) {
    Fidx_base0 <- (cumFindx[idx]) : (cumFindx[idx + 1] - 1)
    FriendIDs <- base0_IDrownums[Fidx_base0 + 1] + 1
    # remove the diag term (idx,idx) from FriendIDs (always the last entry),
    # since we only need to sum all fvec_i[i]^2 once (outside this fun)
    FriendIDs <- FriendIDs[-length(FriendIDs)]
    Dstar_crossprod <- Dstar_crossprod + sum(fvec_i[idx] * fvec_i[FriendIDs])
  }
  return(Dstar_crossprod)
}

# Sum the cross prod vector over connectivity mtx (prod will only appear if (i,j) entry is 1):
est.sigma_sparse <- function(fvec_i, sparse_connectmtx)  {
  n <- length(fvec_i)
  # sum of fvec_i[i]*fvec[j] for correlated cross-product terms (i,j), s.t., i<j
  Dstar_crossprod <- sum_crossprod_Fij(sparseAdjMat = sparse_connectmtx, fvec_i = fvec_i)
  # double cross prod sum + sum of squares over i=1,...,n
  Dstar <- (1/n) * ((2*Dstar_crossprod) + sum(fvec_i^2))
  return(Dstar)
}

est_sigmas <- function(estnames, n, NetInd_k, nF, obsYvals, ests_mat, QY_mat, wts_mat, fWi_mat) {
  `%+%` <- function(a, b) paste0(a, b)

  fWi <- fWi_mat[, "fWi_Qinit"]
  QY.init <- QY_mat[, "QY.init"] 
  h_wts <- wts_mat[, "h_wts"]

  # NetInd as sparse adjacency matrix (new version returns pattern sparse mat ngCMatrix):
  sparse_mat <- NetInd.to.sparseAdjMat(NetInd_k, nF = nF, add_diag = TRUE)
  # Second pass over columns of connectivity mtx to connect indirect intersections (i and j have a common friend but are not friends themselves):
  connectmtx_1stO <- Matrix::crossprod(sparse_mat) # t(sparse_mat)%*%sparse_mat returns nsCMatrix (only non-zero entries)

  # TMLE inference based on the iid IC:
  IC_tmle <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"TMLE",])
  var_tmle <- est.sigma_sparse(IC_tmle, connectmtx_1stO)
  print("var_tmle: " %+% as.numeric(var_tmle/n))

  # ------------------------------------------------------------------------------------------------------------
  # Alternative TMLE variance estimator based on conditional independence of Q(A_i,W_i_ and decomposition of the EIC: 
  # ------------------------------------------------------------------------------------------------------------
  IC_Q_tmle <- h_wts * (obsYvals - QY.init)
  IC_W_tmle <- (fWi - ests_mat[rownames(ests_mat)%in%"TMLE",])
  var_IC_Q_tmle <- (1/n) * sum(IC_Q_tmle^2)
  var_IC_W_tmle <- est.sigma_sparse(IC_W_tmle, connectmtx_1stO)
  var_tmle_2 <- var_IC_Q_tmle + var_IC_W_tmle
  print("var_IC_Q_tmle: " %+% as.numeric(var_IC_Q_tmle/n))
  print("var_IC_W_tmle: " %+% as.numeric(var_IC_W_tmle/n))
  print("var_IC_Q_tmle + var_IC_W_tmle: " %+% as.numeric(var_tmle_2/n))
  # print("total n of non-zero entries in connectmtx_1stO / N^2: "); print(sum(connectmtx_1stO)/(n^2))
  # setting var est. for the data-adaptive target param (conditional on W)
  # var_tmle <- var_IC_Q_tmle

  # ------------------------------------------------------------------------------------------------------------
  # Simple estimator of the iid asymptotic IC-based variance (no adjustment made when two observations i!=j are dependent):
  # ------------------------------------------------------------------------------------------------------------
  iid_var_tmle <- (1/n) * sum(IC_tmle^2)
  # print("iid_var_tmle: " %+% as.numeric(iid_var_tmle/n))

# # ------------------------------------------------------------------------------------------------------------
#   # Alternative estimate, removes everyone who has >25 friends from the connectivity matrix:
#   NetInd_k_2 <- NetInd_k
#   NetInd_k_2[nF >= 20,] <- NA
#   # Ind_Mat <- !is.na(NetInd_k_2)
#   nF_2 <- rowSums(!is.na(NetInd_k_2))
#   sparse_mat_2 <- NetInd.to.sparseAdjMat(NetInd_k_2, nF = nF_2, add_diag = TRUE)
#   connectmtx_1stO_2 <- Matrix::crossprod(sparse_mat_2) # t(sparse_mat)%*%sparse_mat returns nsCMatrix (only non-zero entries)
#   var_tmle <- est.sigma_sparse(IC_tmle, connectmtx_1stO_2)
#   print("total n of non-zero entries in connectmtx_1stO_2 / N^2: "); print(sum(connectmtx_1stO_2)/(n^2))
#   print("var est that excludes all nF >=20"); print(var_tmle/n)
#   # var_tmle_alt <- est.sigma_sparse(IC_tmle, sparse_mat)
#   # var_tmle_alt_2 <- est.sigma_sparse(IC_tmle, sparse_mat_2)
# # ------------------------------------------------------------------------------------------------------------

  # IPTW h (based on the mixture density clever covariate (h)):
  IC_iptw_h <- h_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"h_IPTW",])
  var_iptw_h <- est.sigma_sparse(IC_iptw_h, connectmtx_1stO)
  iid_var_iptw_h <- mean((IC_iptw_h)^2)

  var.ests <- c(abs(var_tmle), abs(var_iptw_h), NA)
  as.var_mat <- matrix(0, nrow = length(var.ests), ncol = 1)
  as.var_mat[,1] <- var.ests
  rownames(as.var_mat) <- estnames; colnames(as.var_mat) <- "Var"

  iid.var.ests = c(iid_var_tmle = abs(iid_var_tmle), # no adjustment for correlations i,j for tmle
                   iid_var_iptw_h = abs(iid_var_iptw_h), # no adjustment for correlations i,j for iptw
                   iid_var_gcomp = NA) 

  iid.vars_mat <- matrix(0, nrow = length(iid.var.ests), ncol = 1)
  iid.vars_mat[,1] <- iid.var.ests
  rownames(iid.vars_mat) <- names(iid.var.ests); colnames(iid.vars_mat) <- "Var"

  # QY.star <- QY_mat[, "QY.star"]
  # g_wts <- wts_mat[,"g_wts"]  
  # var_tmle_A <- var_tmleiptw_1stO <- var_tmleiptw_2ndO <- var_iptw_1stO <- var_iptw_2ndO <- 0
  # var_tmle_A_Q.init <- var_tmle_B_Q.init <- 0
  # # TMLE A (clever covariate update): Inference based on the iid IC analogy, QY.init := initial Q model predictions, h_wts := h_tilde
  # if (!onlyTMLE_B) {
  #   IC_tmle_A <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",])
  #   var_tmle_A <- est.sigma_sparse(IC_tmle_A, connectmtx_1stO)
  # }
  # # TMLE B (weighted model update): Inference based on the iid IC:
  # IC_tmle_B <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",])
  # var_tmle_B <- est.sigma_sparse(IC_tmle_B, connectmtx_1stO)
  # # simple iid estimator of the asymptotic variance (no adjustment made when two observations i!=j are dependent):
  # iid_var_tmle_B <- mean((IC_tmle_B)^2)
  # TMLE based on iptw clever covariate (more non-parametric):
  # if (!onlyTMLE_B) {
  #   IC_tmleiptw <- g_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_g_iptw",])
  #   var_tmleiptw_1stO <- est.sigma_sparse(IC_tmleiptw, connectmtx_1stO)
  # }
  # # IPTW g:
  # if (!onlyTMLE_B) {
  #   iidIC_iptw_g <- g_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"g_iptw",])
  #   var_iptw_1stO <- est.sigma_sparse(iidIC_iptw_g, connectmtx_1stO)
  # }
  # Inference based on the EIC, with factorization into orthogonal components sigma2_DY and sigma2_W_N
  # sigma2_DY_i are independent (since they are conditioned on W,A)
  # sigma2_W_N_i are dependent => need to take double sum of their crossprod among dependent units
  # if (!onlyTMLE_B) {
  #   D_star_Yi.Qinit <- h_wts * (obsYvals - QY.init) # h*(Y-Q_bar_N):
  #   sigma2_DY <- (1/n) * sum(D_star_Yi.Qinit^2)  # Sum_{i} (D_star_Yi)^2

  #   fW_A_i <- fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",]
  #   sigma2_W_N_A <- est.sigma_sparse(fW_A_i, connectmtx_1stO)
  #   var_tmle_A_Q.init <- sigma2_W_N_A + sigma2_DY

  #   # **NEW** TMLE B (weights model update)
  #   fW_B_i <- fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",]
  #   sigma2_W_N_B <- est.sigma_sparse(fW_B_i, connectmtx_1stO)
  #   var_tmle_B_Q.init <- sigma2_W_N_B + sigma2_DY

  #   # D_star_Yi.Qstar <- h_wts * (obsYvals - QY.star)
  #   # D_star_Yi.Qstar[determ.Q] <- 0
  #   # fDY_crossprod <- get.crossprodmtx(D_star_Yi.Qstar)
  #   # double sum over dependent subjects, Sum_{i,j} R_W(i,j)*D_star_Yi*D_star_Yj
  #   # sigma2_Y_N <- est.sigma_sparse(D_star_Yi.Qstar, connectmtx_1stO)
  #   # var_tmle_Q.init_c <- sigma2_W_N_A + sigma2_Y_N

  #   # #--------
  #   # # conservative estimate of the as. variance from EIC for TMLE A:
  #   # # abs terms double sum over dependent subjects, Sum_{i,j} R_W(i,j)*|D_star_Yi|*|D_star_Yj|:
  #   # abs_sigma2_Y_N <- est.sigma_sparse(abs(D_star_Yi.Qstar), connectmtx_1stO)
  #   # var_tmle_A_Q.star_cons <- sigma2_W_N_A + abs_sigma2_Y_N
  #   # # --------
  # }

  # var.ests <- c(abs(var_tmle_A), abs(var_tmle_B), abs(var_tmleiptw_1stO), abs(var_iptw_h), abs(var_iptw_1stO), 0)
  # estnames <- c( "TMLE_A", "TMLE_B", "TMLE_g_IPTW", "h_IPTW", "g_IPTW", "MLE")
  # other.vars = c(
  #               iid_var_tmle_B = abs(iid_var_tmle_B), # no adjustment for correlations i,j
  #               var_tmleiptw_2ndO = abs(var_tmleiptw_2ndO), # adjusting for 2nd order dependence of i,j
  #               var_iptw_2ndO = abs(var_iptw_2ndO), # adjusting for 2nd order dependence of i,j
  #               var_tmle_A_Q.init = abs(var_tmle_A_Q.init), # using the EIC & Q.init for TMLE A
  #               var_tmle_B_Q.init = abs(var_tmle_B_Q.init)  # using the EIC & Q.init for TMLE B
  #               )

  # return(list(as.var_mat = as.var_mat))
  return(list(as.var_mat = as.var_mat, iid.vars_mat = iid.vars_mat))
}


# bootstrap tmle by resampling (sW,sA,Y) with replacement
bootstrap_tmle <- function(estnames, DatNet.ObsP0, tmle_g_out, QY_mat, wts_mat, fWi_mat) {
  # browser()
  # names(tmle_g_out)

  # QY.init <- DatNet.ObsP0$noNA.Ynodevals # getting all node vals, inc. deterministic  
  # QY.init[!DatNet.ObsP0$det.Y] <- m.Q.init$predict(newdata = DatNet.ObsP0)$getprobA1[!DatNet.ObsP0$det.Y] # getting predictions P(Y=1) for non-DET Y
  QY.init <- QY_mat[, "QY.init"] 
  off <- qlogis(QY.init)  # offset

  # fWi <- fWi_mat[, "fWi_Qinit"]

  # fit.hbars_t <- system.time(fit.hbars.res <- fit.hbars(DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_params_list)) # fit the clever covariate
  # m.h.fit <- fit.hbars.res$m.h.fit
  DatNet.gstar <- tmle_g_out$DatNet.gstar

  psi.evaluator <- tmle_g_out$psi.evaluator
  m.Q.init <- tmle_g_out$m.Q.init
  m.h.fit <- tmle_g_out$m.h.fit
  
  # h_wts <- tmle_g_out$h_wts
  h_wts <- wts_mat[, "h_wts"]

  Y <- DatNet.ObsP0$noNA.Ynodevals
  determ.Q <- DatNet.ObsP0$det.Y

  # m.h.fit$summeas.g0,
  # m.h.fit$summeas.gstar

  #************************************************
  # BOOTSTRAP TMLE UPDATES:
  #************************************************
  n_boots <- 100
  boot_eps <- vector(mode = "numeric", length = n_boots)
  boot_tmle_B <- vector(mode = "numeric", length = n_boots)
  boot_time <- system.time(
    for (i in (1:n_boots)) {
        boot_idx <- sample.int(n = DatNet.ObsP0$nobs, replace = TRUE)
        boot.tmle.obj <- tmle.update(estnames = estnames,
                                     Y = Y[boot_idx], off = off[boot_idx], h_wts = h_wts[boot_idx], 
                                     determ.Q = determ.Q[boot_idx], predictQ = FALSE)
        boot_eps[i] <- boot.tmle.obj$m.Q.star.coef
        boot_tmle_B[i] <- mean(psi.evaluator$get.boot.tmleB(m.Q.starB = boot_eps[i], boot_idx = boot_idx))
    }
  )
  var_tmleB_boot <- var(boot_tmle_B)

  print("boot_time for n_boots = " %+% n_boots); print(boot_time)
  print("mean(boot_eps): "); print(mean(boot_eps))
  print("mean(boot_tmle_B): "); print(mean(boot_tmle_B))
  print("var(boot_tmle_B): "); print(var(boot_tmle_B))

  return(var_tmleB_boot)

}

# parametric bootstrap
# f.g0, 
par_bootstrap_tmle <- function(estnames, DatNet.ObsP0, tmle_g_out, QY_mat, wts_mat, fWi_mat) {
  f.g0 <- function(data) {
    return(A = rbinom(n = nrow(data), size = 1, prob = 0.25))
  }

  # names(tmle_g_out)

  n_boots <- 100
  boot_eps <- vector(mode = "numeric", length = n_boots)
  boot_gcomp <- vector(mode = "numeric", length = n_boots)
  boot_tmle_B <- vector(mode = "numeric", length = n_boots)

  # 0. Save the original input data.table OdataDT, otherwise it will be over-written

  psi.evaluator <- tmle_g_out$psi.evaluator
  
  if (!DatNet.ObsP0$datnetW$Odata$curr.data.A.g0) DatNet.ObsP0$datnetW$Odata$restoreAnodes() # use the observed Anodes

  OdataDT <- DatNet.ObsP0$datnetW$Odata$OdataDT

  DatNet.g0.boot <- DatNet.ObsP0
  DatNet.gstar <- tmle_g_out$DatNet.gstar

# OdataDT
#           ID HUB PA nF nF.PA W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 A sum.net.A        probY Y
#     1:     1   1  0 46     4  4  1      0        1        0        0        1        1 1         6 0.9992165879 1
#     2:     2   1  0 56     8  2  0      1        1        1        1        1        1 1         3 0.9517114768 1
#     3:     3   1  0 49     8  4  0      0        0        1        0        1        0 1         5 0.9998742474 1
#     4:     4   1  0 47     6  5  0      0        0        0        0        0        0 1         1 0.1456751954 0
#     5:     5   1  1 54     5  4  0      0        1        1        1        1        0 1         5 0.9999999927 1
#    ---                                                                                                           
# 49996: 49996   0  0  5     0  4  0      0        0        0        1        1        0 0         0 0.0005424175 0
# 49997: 49997   0  0  5     0  4  1      1        1        0        1        1        0 0         2 0.0431513304 1
# 49998: 49998   0  0  5     2  3  0      0        1        1        1        0        1 0         0 0.0001210807 0
# 49999: 49999   0  0  5     0  2  0      0        0        0        1        0        0 0         3 0.0001554654 0
# 50000: 50000   0  0  5     0  3  0      1        1        1        1        0        1 0         2 0.0001554654 0

# mean(OdataDT[["A"]]) # A is under g.star now
# [1] 0.09728

# head(DatNet.g0.boot$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,]   1   1  4  1      0        1        0        0        1        1 46 1     4       1             0         6                  24
# [2,]   1   1  2  0      1        1        1        1        1        1 56 0     8       0             0         3                  24
# [3,]   1   1  4  0      0        0        1        0        1        0 49 0     8       0             0         5                  40
# [4,]   1   1  5  0      0        0        0        0        0        0 47 1     6       1             0         1                   6
# [5,]   1   0  4  0      0        1        1        1        1        0 54 0     5       0             0         5                  25
# [6,]   1   1  4  0      1        0        1        0        0        1 54 0     5       0             0         3                  15

# head(DatNet.gstar$dat.sVar)
#     HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,]   1   1  4  1      0        1        0        0        1        1 46 1     4       1             0        11                  44
# [2,]   1   1  2  0      1        1        1        1        1        1 56 1     8       1             0        13                 104
# [3,]   1   1  4  0      0        0        1        0        1        0 49 1     8       1             0        12                  96
# [4,]   1   1  5  0      0        0        0        0        0        0 47 1     6       1             0         4                  24
# [5,]   1   0  4  0      0        1        1        1        1        0 54 1     5       0             0         5                  25
# [6,]   1   1  4  0      1        0        1        0        0        1 54 1     5       1             0         7                  35

# head(DatNet.g0.boot$datnetW$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF
# [1,]   1   1  4  1      0        1        0        0        1        1 46
# [2,]   1   1  2  0      1        1        1        1        1        1 56
# [3,]   1   1  4  0      0        0        1        0        1        0 49
# [4,]   1   1  5  0      0        0        0        0        0        0 47
# [5,]   1   0  4  0      0        1        1        1        1        0 54
# [6,]   1   1  4  0      1        0        1        0        0        1 54

# head(DatNet.gstar$datnetW$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF
# [1,]   1   1  4  1      0        1        0        0        1        1 46
# [2,]   1   1  2  0      1        1        1        1        1        1 56
# [3,]   1   1  4  0      0        0        1        0        1        0 49
# [4,]   1   1  5  0      0        0        0        0        0        0 47
# [5,]   1   0  4  0      0        1        1        1        1        0 54
# [6,]   1   1  4  0      1        0        1        0        0        1 54

# head(DatNet.g0.boot$datnetA$dat.sVar)
#      A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,] 1     4       1             0         6                  24
# [2,] 0     8       0             0         3                  24
# [3,] 0     8       0             0         5                  40
# [4,] 1     6       1             0         1                   6
# [5,] 0     5       0             0         5                  25
# [6,] 0     5       0             0         3                  15

# head(DatNet.gstar$datnetA$dat.sVar)
#      A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,] 1     4       1             0         6                  24
# [2,] 0     8       0             0         3                  24
# [3,] 0     8       0             0         5                  40
# [4,] 1     6       1             0         1                   6
# [5,] 0     5       0             0         5                  25
# [6,] 0     5       0             0         3                  15


# loop over n_boots
  boot_time <- system.time(
    for (i in (1:n_boots)) {

  # 1. Resample W (with replacement) by re-purposing the instance of DatNet.gstar. Re-shuffle pre-saved values of Y and det.Y
    boot_idx <- sample.int(n = DatNet.ObsP0$nobs, replace = TRUE)
    DatNet.g0.boot$datnetW$Odata$OdataDT <- OdataDT[boot_idx, ]
    DatNet.g0.boot$datnetW$make.sVar(sVar.object = tmle_g_out$sW)

    DatNet.g0.boot$det.Y <- DatNet.g0.boot$det.Y[boot_idx]
    DatNet.g0.boot$noNA.Ynodevals <- DatNet.g0.boot$noNA.Ynodevals[boot_idx]

# head(DatNet.g0.boot$datnetW$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF
# [1,]   0   1  2  0      0        1        1        1        1        1 46
# [2,]   0   1  4  1      1        0        0        0        1        0 56
# [3,]   0   1  3  0      1        0        1        0        1        1 49
# [4,]   0   1  1  0      1        1        1        1        0        1 47
# [5,]   0   1  3  1      1        1        1        0        0        0 54
# [6,]   0   1  4  1      0        1        0        1        1        1 54
# head(DatNet.gstar$datnetW$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF
# [1,]   0   1  2  0      0        1        1        1        1        1 46
# [2,]   0   1  4  1      1        0        0        0        1        0 56
# [3,]   0   1  3  0      1        0        1        0        1        1 49
# [4,]   0   1  1  0      1        1        1        1        0        1 47
# [5,]   0   1  3  1      1        1        1        0        0        0 54
# [6,]   0   1  4  1      0        1        0        1        1        1 54


# head(DatNet.g0.boot$datnetW$Odata$OdataDT)
#       ID HUB PA nF nF.PA W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 A sum.net.A        probY Y
# 1: 26690   0  0  6     0  2  0      0        1        1        1        1        1 0         4 1.473071e-03 0
# 2: 37923   0  0  7     1  4  1      1        0        0        0        1        0 0         3 4.315133e-02 0
# 3:  3414   0  0 23     2  3  0      1        0        1        0        1        1 0         3 2.426364e-03 0
# 4: 27356   0  0 11     1  1  0      1        1        1        1        0        1 0         1 5.719814e-05 0
# 5: 49128   0  0  5     0  3  1      1        1        1        0        0        0 0         3 6.066238e-03 0
# 6: 16666   0  0 13     1  4  1      0        1        0        1        1        1 0         2 2.085815e-02 1
# head(DatNet.gstar$datnetW$Odata$OdataDT)
#       ID HUB PA nF nF.PA W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 A sum.net.A        probY Y
# 1: 26690   0  0  6     0  2  0      0        1        1        1        1        1 0         4 1.473071e-03 0
# 2: 37923   0  0  7     1  4  1      1        0        0        0        1        0 0         3 4.315133e-02 0
# 3:  3414   0  0 23     2  3  0      1        0        1        0        1        1 0         3 2.426364e-03 0
# 4: 27356   0  0 11     1  1  0      1        1        1        1        0        1 0         1 5.719814e-05 0
# 5: 49128   0  0  5     0  3  1      1        1        1        0        0        0 0         3 6.066238e-03 0
# 6: 16666   0  0 13     1  4  1      0        1        0        1        1        1 0         2 2.085815e-02 1
# head(DatNet.gstar$datnetW$Odata$OdataDT)
#       ID HUB PA nF nF.PA W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 A sum.net.A        probY Y
# 1: 26690   0  0  6     0  2  0      0        1        1        1        1        1 1         4 1.473071e-03 0
# 2: 37923   0  0  7     1  4  1      1        0        0        0        1        0 1         3 4.315133e-02 0
# 3:  3414   0  0 23     2  3  0      1        0        1        0        1        1 1         3 2.426364e-03 0
# 4: 27356   0  0 11     1  1  0      1        1        1        1        0        1 1         1 5.719814e-05 0
# 5: 49128   0  0  5     0  3  1      1        1        1        0        0        0 1         3 6.066238e-03 0
# 6: 16666   0  0 13     1  4  1      0        1        0        1        1        1 1         2 2.085815e-02 1
  
  # 2. Generate new A's from g0 or g.N (replace A with sampled A's in DatNet.ObsP0) & Re-create the summary measures (sW,sA) based on new DatNet.ObsP0

    # head(DatNet.g0.boot$datnetA$dat.sVar)
    # datnetA_make.sVar_time <-  system.time(
    #   DatNet.g0.boot$datnetA$make.sVar(sVar.object = tmle_g_out$sA)
    #   )
    # print("datnetA_make.sVar_time"); print(datnetA_make.sVar_time)    
    # head(DatNet.g0.boot$datnetA$dat.sVar)

    # head(DatNet.g0.boot$dat.sVar)
    # make.dat.sWsA_time <- system.time(
      DatNet.g0.boot$make.dat.sWsA(p = 1, f.g_fun = f.g0, sA.object = tmle_g_out$sA, DatNet.ObsP0 = DatNet.g0.boot)
      # )
    # print("make.dat.sWsA_time: "); print(make.dat.sWsA_time)

# head(DatNet.g0.boot$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,]   0   1  2  0      0        1        1        1        1        1 46 0     0       0             0        10                   0
# [2,]   0   1  4  1      1        0        0        0        1        0 56 0     1       0             0        18                  18
# [3,]   0   1  3  0      1        0        1        0        1        1 49 0     2       0             0        17                  34
# [4,]   0   1  1  0      1        1        1        1        0        1 47 1     1       1             0         9                   9
# [5,]   0   1  3  1      1        1        1        0        0        0 54 1     0       1             0        12                   0
# [6,]   0   1  4  1      0        1        0        1        1        1 54 0     1       0             0        15                  15
# head(DatNet.gstar$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,]   1   1  4  1      0        1        0        0        1        1 46 1     4       1             0        11                  44
# [2,]   1   1  2  0      1        1        1        1        1        1 56 1     8       1             0        13                 104
# [3,]   1   1  4  0      0        0        1        0        1        0 49 1     8       1             0        12                  96
# [4,]   1   1  5  0      0        0        0        0        0        0 47 1     6       1             0         4                  24
# [5,]   1   0  4  0      0        1        1        1        1        0 54 1     5       0             0         5                  25
# [6,]   1   1  4  0      1        0        1        0        0        1 54 1     5       1             0         7                  35
# after re-creating DatNet.gstar under bootstrapped W:
# head(DatNet.gstar$dat.sVar)
#      HUB PA0 W1 W2 WNoise corrW.F1 corrW.F2 corrW.F3 corrW.F4 corrW.F5 nF A nF.PA A.PAeq0 nFPAeq0.PAeq1 sum.net.A sum.net.A.sum.netPA
# [1,]   0   1  2  0      0        1        1        1        1        1 46 1     0       1             0        36                   0
# [2,]   0   1  4  1      1        0        0        0        1        0 56 1     1       1             0        37                  37
# [3,]   0   1  3  0      1        0        1        0        1        1 49 1     2       1             0        37                  74
# [4,]   0   1  1  0      1        1        1        1        0        1 47 1     1       1             0        37                  37
# [5,]   0   1  3  1      1        1        1        0        0        0 54 1     0       1             0        33                   0
# [6,]   0   1  4  1      0        1        0        1        1        1 54 1     1       1             0        32                  32

    # DatNet.g0 <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
    # head(DatNet.g0.boot$dat.sVar)
    # DatNet.g0.boot$make.dat.sWsA()
    # head(DatNet.g0.boot$dat.sVar)
    # mean(DatNet.g0.boot$dat.sVar[,"A"])
    # DatNet.g0.boot$datnetW$Odata$OdataDT
    # mean(DatNet.g0.boot$datnetW$Odata$OdataDT[["A"]])
    # #---------------------------------------------------------------------------------
    # # BUILDING OBSERVED sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
    # #---------------------------------------------------------------------------------
    # datnetW <- DatNet$new(netind_cl = netind_cl, nFnode = nFnode)
    # # datnetW$make.sVar(Odata = data, sVar.object = sW)
    # datnetW$make.sVar(Odata = OdataDT_R6[], sVar.object = sW)
    # datnetW$fixmiss_sVar() # permanently replace NA values in sW with 0
    # datnetA <- DatNet$new(netind_cl = netind_cl)
    # # datnetA$make.sVar(Odata = data, sVar.object = sA)
    # datnetA$make.sVar(Odata = OdataDT_R6, sVar.object = sA)
    # DatNet.ObsP0 <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()

  # 4. Predict P(Y_i=1|sW,sA) using m.Q.init (the initial fit \bar{Q}_N) based on newly resampled (sW,sA)

    # these are the new bootstrapped Y's that are deterministically assigned:
    detY.boot <- DatNet.ObsP0$det.Y
    QY.init.boot <- DatNet.ObsP0$noNA.Ynodevals
    QY.init.boot[!detY.boot] <- tmle_g_out$m.Q.init$predict(newdata = DatNet.g0.boot)$getprobA1[!detY.boot] # getting predictions P(Y=1) for non-DET Y
    off.boot <- qlogis(QY.init.boot)  # offset
    # length(QY.init.boot); head(QY.init.boot)

  # 5. Sample a vector of new (Y_i, i=1,...,N)

    # psi.evaluator$sampleY(Qprob = QY.init)?
    Y.boot <- rbinom(n = length(QY.init.boot), size = 1, prob = QY.init.boot)
    # length(Y.boot); head(Y.boot)

  # 6. Predict new weights h_wts = P_{\bar{g}^*_N}(sA | sW)/P_{\bar{g}_0}(sA | sW) 
  # using the existing fits m.h.fit$summeas.g0 and m.h.fit$summeas.gstar and the newly resampled (sW,sA)

    h_wts.boot <- predict.hbars(newdatnet = DatNet.g0.boot, m.h.fit = tmle_g_out$m.h.fit)
    # head(h_wts.boot)

  # 7. Fit a TMLE update epsilon on this new data.

    boot.tmle.obj <- tmle.update(estnames = estnames,
                                 Y = Y.boot, off = off.boot, h_wts = h_wts.boot, 
                                 determ.Q = detY.boot, predictQ = FALSE)

    boot_eps[i] <- boot.tmle.obj$m.Q.star.coef


  # 8. Re-create DatNet.gstar with boostrapped summaires sW and sA generated under f.gstar

    DatNet.gstar$make.dat.sWsA(p = 1, f.g_fun = tmle_g_out$f.gstar, sA.object = tmle_g_out$sA, DatNet.ObsP0 = DatNet.g0.boot)

  # 9. Evaluate the substitution estimator and the components of the EIC D_Y and D_W
    # head(DatNet.gstar$dat.sVar)
    # ******** STILL NEED TO VERIFY THAT get.gcomp & get.tmleB IS ACTUALLY USING BOOTSTRAPPED DATA....*********
    boot_gcomp[i] <- mean(psi.evaluator$get.gcomp(m.Q.init = tmle_g_out$m.Q.init))
    boot_tmle_B[i] <- mean(psi.evaluator$get.tmleB(m.Q.starB = boot.tmle.obj$m.Q.star.coef))

    }
  )

  # browser()
  var_gcomp_boot <- var(boot_gcomp)
  var_tmleB_boot <- var(boot_tmle_B)
  
  # plot(density(boot_gcomp))
  # plot(density(boot_tmle_B))

  print("boot_time for n_boots = " %+% n_boots); print(boot_time)
  #  [1] "boot_time for n_boots = 10"
  #   user  system elapsed 
  # 35.144  12.216  51.154 
  print("mean(boot_eps): "); print(mean(boot_eps))
  print("mean(boot_tmle_B): "); print(mean(boot_tmle_B))
  print("var(boot_tmle_B): "); print(var(boot_tmle_B))
  print("var(var_gcomp_boot): "); print(var(boot_gcomp))

  return(var_tmleB_boot)
}


# create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
make_EYg_obj <- function(estnames, estoutnames, alpha, DatNet.ObsP0, tmle_g_out, tmle_g2_out=NULL) {
  nobs <- DatNet.ObsP0$nobs
  NetInd_k <- DatNet.ObsP0$netind_cl$NetInd_k
  nF <- DatNet.ObsP0$netind_cl$nF
  ests_mat <- tmle_g_out$ests_mat
  QY_mat <- tmle_g_out$QY_mat
  fWi_mat <- tmle_g_out$fWi_mat
  wts_mat <- tmle_g_out$wts_mat

  if (!is.null(tmle_g2_out)) {
    ests_mat <- tmle_g_out$ests_mat - tmle_g2_out$ests_mat
    fWi_mat <- tmle_g_out$fWi_mat - tmle_g2_out$fWi_mat
    wts_mat <- tmle_g_out$wts_mat - tmle_g2_out$wts_mat
  }

  # ------------------------------------------------------------------------------------------
  # get the IC-based asymptotic variance estimates:
  # ------------------------------------------------------------------------------------------
  getVar_time <- system.time(
    as.vars_obj <- est_sigmas(estnames = estnames, n = nobs, 
                              NetInd_k = NetInd_k, nF = nF,
                              obsYvals = DatNet.ObsP0$noNA.Ynodevals,
                              ests_mat = ests_mat, QY_mat = QY_mat, 
                              wts_mat = wts_mat, fWi_mat = fWi_mat)
  )
  if (gvars$verbose) {
    print("time to estimate Vars: "); print(getVar_time)  
  }

  get_CI <- function(xrow, n) {
    f_est_CI <- function(n, psi, sigma2_N) { # get CI
      z_alpha <- qnorm(1-alpha/2)
      CI_est <- c(psi - z_alpha*sqrt(sigma2_N) / sqrt(n), psi + z_alpha*sqrt(sigma2_N) / sqrt(n))
      return(CI_est)
    }
    psi <- xrow["estimate"];
    sigma2_N <- xrow["Var"];
    return(f_est_CI(n = n, psi = psi, sigma2_N = sigma2_N))
  }

  CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$as.var_mat), 1, get_CI, n = nobs))
  colnames(CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))
  
  # CIs based on IID variance :
  iid.CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$iid.vars_mat), 1, get_CI, n = nobs))
  colnames(iid.CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))


  # ------------------------------------------------------------------------------------------
  # PARAMETRIC BOOSTRAP:
  # ------------------------------------------------------------------------------------------
  # var_tmleB_boot <- bootstrap_tmle(estnames, DatNet.ObsP0, tmle_g_out, QY_mat, wts_mat, fWi_mat)
  var_tmleB_boot <- par_bootstrap_tmle(estnames, DatNet.ObsP0, tmle_g_out, QY_mat, wts_mat, fWi_mat)


  # ------------------------------------------------------------------------------------------
  # RENAME ESTIMATORS FOR THE FINAL OUTPUT:
  # ------------------------------------------------------------------------------------------
  rownames(ests_mat) <- estoutnames
  rownames(as.vars_obj$as.var_mat) <- estoutnames
  rownames(as.vars_obj$iid.vars_mat) <- estoutnames
  rownames(CIs_mat) <- estoutnames
  rownames(iid.CIs_mat) <- estoutnames

  EY_g.star <- list(estimates = ests_mat,
                    vars = (as.vars_obj$as.var_mat / nobs),
                    CIs = CIs_mat,
                    iid.vars = (as.vars_obj$iid.vars_mat / nobs),
                    # var_tmleB_boot = tmle_g_out$var_tmleB_boot,
                    var_tmleB_boot = var_tmleB_boot,
                    iid.CIs = iid.CIs_mat,
                    h_g0_SummariesModel = NULL,
                    h_gstar_SummariesModel = NULL)

  if (is.null(tmle_g2_out)) {
    EY_g.star[["h_g0_SummariesModel"]] <- tmle_g_out$h_g0_SummariesModel
    EY_g.star[["h_gstar_SummariesModel"]] <- tmle_g_out$h_gstar_SummariesModel
  }

  return(EY_g.star)
}