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
sum_crossprod_Fij <- function(sparseAdjMat, fvec_i, excludeIDs = NULL) {
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
    # always remove the diag term (idx,idx) from FriendIDs (always the last entry),
    # since we only need to sum all fvec_i[i]^2 once (outside this fun)
    if (!is.null(excludeIDs)) {
      exclude_now <- c(which(FriendIDs %in% excludeIDs), length(FriendIDs))
    } else {
      exclude_now <- length(FriendIDs)
    }
    # FriendIDs <- FriendIDs[-length(FriendIDs)]
    FriendIDs <- FriendIDs[-exclude_now]
    Dstar_crossprod <- Dstar_crossprod + sum(fvec_i[idx] * fvec_i[FriendIDs])
  }
  return(Dstar_crossprod)
}

# Sum the cross prod vector over connectivity mtx (prod will only appear if (i,j) entry is 1):
est.sigma_sparse <- function(fvec_i, sparse_connectmtx, excludeIDs = NULL)  {
  n <- length(fvec_i)
  # sum of fvec_i[i]*fvec[j] for correlated cross-product terms (i,j), s.t., i<j
  Dstar_crossprod <- sum_crossprod_Fij(sparseAdjMat = sparse_connectmtx, fvec_i = fvec_i, excludeIDs = excludeIDs)
  # double cross prod sum + sum of squares over i=1,...,n
  Dstar <- (1/n) * ((2*Dstar_crossprod) + sum(fvec_i^2))
  return(Dstar)
}

est_sigmas <- function(estnames, n, NetInd_k, nF, obsYvals, ests_mat, QY_mat, wts_mat, fWi_mat, MC.tmle.eval) {
  `%+%` <- function(a, b) paste0(a, b)
  fWi <- fWi_mat[, "fWi_Qinit"]
  QY.init <- QY_mat[, "QY.init"]
  h_wts <- wts_mat[, "h_wts"]

  # NetInd as sparse adjacency matrix (new version returns pattern sparse mat ngCMatrix):
  sparse_mat <- NetInd.to.sparseAdjMat(NetInd_k, nF = nF, add_diag = TRUE)
  # Second pass over columns of connectivity mtx to connect indirect intersections (i and j have a common friend but are not friends themselves):
  connectmtx_1stO <- Matrix::crossprod(sparse_mat) # t(sparse_mat)%*%sparse_mat returns nsCMatrix (only non-zero entries)

  # observation-specific IC values:
  IC_tmle <- h_wts * (obsYvals - QY.init) + (fWi - MC.tmle.eval)
  # IC_tmle <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"TMLE",])

  # TMLE variance estimator based on the above IC that adjusts for correlated observation (a double sum):
  var_tmle <- est.sigma_sparse(IC_tmle, connectmtx_1stO)

  # ------------------------------------------------------------------------------------------------------------
  # Alternative TMLE variance estimators based on conditional independence of Q(A_i,W_i_ and decomposition of the EIC:
  # var_IC_Q_tmle gives inference CONDITIONAL all W.
  # ------------------------------------------------------------------------------------------------------------
  n_knockout <- c(20, 200, 2000)
  aux.vars_mat <- matrix(0, nrow = 1+length(n_knockout), ncol = 1)
  colnames(aux.vars_mat) <- "Var"

  IC_Q_tmle <- h_wts * (obsYvals - QY.init)
  var_IC_Q_tmle <- (1/n) * sum(IC_Q_tmle^2)

  IC_W_tmle <- (fWi - MC.tmle.eval)
  # IC_W_tmle <- (fWi - ests_mat[rownames(ests_mat)%in%"TMLE",])

  var_IC_W_tmle <- est.sigma_sparse(IC_W_tmle, connectmtx_1stO)
  var_tmle_2 <- var_IC_Q_tmle + var_IC_W_tmle
  var_tmle_noindirect <- est.sigma_sparse(IC_tmle, sparse_mat)
  aux.vars_mat[1,1] <- var_tmle_2
  for (ID_idx in n_knockout) {
    aux.vars_mat[which(n_knockout %in% ID_idx) + 1, 1] <- est.sigma_sparse(IC_tmle, connectmtx_1stO, excludeIDs = c(1:ID_idx))
  }
  rownames(aux.vars_mat) <- c("factorized_DY_DW", "noID_1_to_" %+% n_knockout)

  # ------------------------------------------------------------------------------------------------------------
  # Simple estimator of the iid asymptotic IC-based variance (no adjustment made when two observations i!=j are dependent):
  # ------------------------------------------------------------------------------------------------------------
  iid_var_tmle <- (1/n) * sum(IC_tmle^2)
  IC_vars <- rbind(var_tmle/n, var_IC_Q_tmle/n, iid_var_tmle/n)
  rownames(IC_vars) <- c("var_IC_tmle", "var_IC_Q_tmle", "iid_var_tmle")
  # print(IC_vars)

  # IPTW h (based on the mixture density clever covariate (h)):
  IC_iptw_h <- h_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"h_IPTW",])
  var_iptw_h <- est.sigma_sparse(IC_iptw_h, connectmtx_1stO)
  iid_var_iptw_h <- mean((IC_iptw_h)^2)

  var.ests <- c(abs(var_tmle), abs(var_iptw_h), NA)
  as.var_mat <- matrix(0, nrow = length(var.ests), ncol = 1)
  as.var_mat[,1] <- var.ests
  rownames(as.var_mat) <- estnames; colnames(as.var_mat) <- "Var"

  # Inference conditional on all W
  condW.vars.ests <- c(condW_var_tmle = abs(var_IC_Q_tmle),
                      condW_var_iptw_h = NA,
                      condW_var_gcomp = NA)

  condW.vars_mat <- matrix(0, nrow = length(condW.vars.ests), ncol = 1)
  condW.vars_mat[,1] <- condW.vars.ests
  rownames(condW.vars_mat) <- names(condW.vars.ests); colnames(condW.vars_mat) <- "Var"

  # print("condW.vars_mat: "); print(condW.vars_mat/n)

  # IID inference (ignores all dependence)
  iid.vars.ests = c(iid_var_tmle = abs(iid_var_tmle), # no adjustment for correlations i,j for tmle
                    iid_var_iptw_h = abs(iid_var_iptw_h), # no adjustment for correlations i,j for iptw
                    iid_var_gcomp = NA)
  iid.vars_mat <- matrix(0, nrow = length(iid.vars.ests), ncol = 1)
  iid.vars_mat[,1] <- iid.vars.ests
  rownames(iid.vars_mat) <- names(iid.vars.ests); colnames(iid.vars_mat) <- "Var"

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
  return(list(as.var_mat = as.var_mat, condW.vars_mat = condW.vars_mat, iid.vars_mat = iid.vars_mat, aux.vars_mat = aux.vars_mat))
}

# recursively check for outvar name until found, then return the appropriate SummaryModel object:
findRegSummaryObj <- function(fit.obj, outvar) {
  if (is.list(fit.obj$outvar) | (length(fit.obj$outvar)>1)) {
    return(unlist(lapply(fit.obj$getPsAsW.models(), findRegSummaryObj, outvar)))
  } else if (fit.obj$outvar %in% outvar) {
    return(fit.obj)
  }
}

fit_reg.forms <- function(boot.form, DatNet.ObsP0, sW, sA) {
  boot.form.sVars <- process_regforms(boot.form, sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
  pred_nms <- boot.form.sVars$predvars
  boot_outvar_nms <- boot.form.sVars$outvars
  check.pred.exist <- unlist(lapply(unlist(pred_nms), function(sWname) sWname %in% c(DatNet.ObsP0$datnetW$names.sVar,DatNet.ObsP0$datnetA$names.sVar)))
  if (!all(check.pred.exist)) stop("the following predictors from boot.form regression could not be located in sW/sA summary measures: " %+%
                                    paste0(pred_nms[!check.pred.exist], collapse = ","))
  check.outvar.exist <- unlist(lapply(boot_outvar_nms, function(sAname) sAname %in% c(DatNet.ObsP0$datnetW$names.sVar,DatNet.ObsP0$datnetA$names.sVar)))
  if (!all(check.outvar.exist)) stop("the following outcomes from boot.form regression could not be located in sW/sA summary measures: " %+%
                                    paste0(boot_outvar_nms[!check.outvar.exist], collapse = ","))
  outvar_class <- lapply(boot_outvar_nms, function(sA_nms) DatNet.ObsP0$datnetA$type.sVar[sA_nms])
  subsets_expr <- lapply(boot_outvar_nms, function(var) lapply(var, function(var) {var}))
  regclass.boot <- RegressionClass$new(sep_predvars_sets = TRUE,
                                     outvar.class = outvar_class,
                                     outvar = boot_outvar_nms,
                                     predvars = pred_nms,
                                     subset = subsets_expr)
  regclass.boot$S3class <- "generic"
  summeas_boot <- newsummarymodel(reg = regclass.boot, DatNet.sWsA.g0 = DatNet.ObsP0)$fit(data = DatNet.ObsP0)
  return(list(boot_outvar_nms = boot_outvar_nms, summeas_boot = summeas_boot))
}


# --------------------------------------------------------------------------------------------------------
# Evaluating the D_Wi component of the EIC IC with Monte-Carlo sampling (requires evaluation of psi_i_n)
# --------------------------------------------------------------------------------------------------------
MCeval_fWi <- function(n.MC, DatNet.ObsP0, tmle_g1_out, tmle_g2_out) {
  mean_obs_gcomp_g1 <- vector(mode = "numeric", length = DatNet.ObsP0$nobs)
  mean_obs_gcomp_g2 <- vector(mode = "numeric", length = DatNet.ObsP0$nobs)
  mean_obs_tmle_B_g1 <- vector(mode = "numeric", length = DatNet.ObsP0$nobs)
  mean_obs_tmle_B_g2 <- vector(mode = "numeric", length = DatNet.ObsP0$nobs)

  mean_obs_gcomp_g1[] <- mean_obs_gcomp_g2[] <- mean_obs_tmle_B_g1[] <- mean_obs_tmle_B_g2[] <- 0L

  psi.evaluator <- tmle_g1_out$psi.evaluator
  # -----------------------------------------------------------------------------------------------
  # Always start with the observed Anodes
  # Save the original input data.table OdataDT, otherwise it will be over-written:
  # Note we are not saving the original saved Anodes and sA -> these fields NULLed during bootstrap
  # -----------------------------------------------------------------------------------------------
  if (!DatNet.ObsP0$Odata$curr_data_A_g0) DatNet.ObsP0$Odata$restoreAnodes()
  OdataDT.P0 <- DatNet.ObsP0$Odata$OdataDT
  noNA.Ynodevals.P0 <- DatNet.ObsP0$noNA.Ynodevals
  det.Y.P0 <- DatNet.ObsP0$det.Y
  DatNet.gstar <- tmle_g1_out$DatNet.gstar

  # -----------------------------------------------------------------------------------------------
  # loop over n.MC
  # -----------------------------------------------------------------------------------------------
  for (i in (1:n.MC)) {
    # 1. Resample W (with replacement) as if IID (re-purposing the instance of DatNet.gstar); Re-evaluate baseline summaries sW;
    boot_idx <- sample.int(n = DatNet.ObsP0$nobs, replace = TRUE)
    DatNet.ObsP0$Odata$OdataDT <- OdataDT.P0[boot_idx, ]
    # Need to NULL previously backed-up values of A and sA, since they no longer correspond with bootstrapped sample
    DatNet.ObsP0$Odata$A_g0_DT <- NULL
    DatNet.ObsP0$Odata$sA_g0_DT <- NULL
    DatNet.ObsP0$Odata$save_sA_Vars <- NULL
    DatNet.ObsP0$datnetW$make.sVar(sVar.object = tmle_g1_out$sW)
    DatNet.ObsP0$datnetW$fixmiss_sVar() # permanently replace NA values in sW with 0
    DatNet.ObsP0$det.Y <- det.Y.P0[boot_idx]
    DatNet.ObsP0$noNA.Ynodevals <- noNA.Ynodevals.P0[boot_idx]

    # 2. USE old observed A's (no re-sampling)
    for (Anode in DatNet.ObsP0$nodes$Anodes) {
      DatNet.ObsP0$Odata$replaceOneAnode(AnodeName = Anode, newAnodeVal = OdataDT.P0[[Anode]])
    }
    # Re-evaluate exposure summaries based on resampled W and old (observed) A's
    DatNet.ObsP0$datnetA$make.sVar(sVar.object = tmle_g1_out$sA)
    DatNet.ObsP0$Odata$curr_data_A_g0 <- TRUE

    # 4. Predict P(Y_i=1|sW,sA) using m.Q.init (the initial fit \bar{Q}_N) based on newly resampled W and re-computed summaries (sW,sA)
    detY.boot <- DatNet.ObsP0$det.Y
    QY.init.boot <- DatNet.ObsP0$noNA.Ynodevals
    QY.init.boot[!detY.boot] <- tmle_g1_out$m.Q.init$predict(newdata = DatNet.ObsP0)$getprobA1[!detY.boot] # getting predictions P(Y=1) for non-DET Y
    # off.boot <- qlogis(QY.init.boot)  # offset

    # 8. Re-create DatNet.gstar with boostrapped W and sA generated under f.gstar:
    # However, first need to save (Anodes, sA) from the observed sample, otherwise they are forever lost
    DatNet.ObsP0$Odata$backupAnodes(sA = tmle_g1_out$sA)

    # Will over-write Anodes/sA in DatNet.ObsP0$Odata with tmle_g1_out$new.sA:
    DatNet.gstar$make.dat.sWsA(p = 1, new.sA.object = tmle_g1_out$new.sA, sA.object = tmle_g1_out$sA, DatNet.ObsP0 = DatNet.ObsP0)

    # 9. Evaluate the substitution estimator and the components of the EIC D_Y and D_W:
    fWi.g1 <- psi.evaluator$get.gcomp(m.Q.init = tmle_g1_out$m.Q.init)
    tmle.g1 <- psi.evaluator$get.tmleB(m.Q.starB = tmle_g1_out$MC_fit_params$m.Q.star)

    mean_obs_gcomp_g1 <- mean_obs_gcomp_g1 + fWi.g1
    mean_obs_tmle_B_g1 <- mean_obs_tmle_B_g1 + tmle.g1

    # 10. If (!is.null(tmle_g2_out)) then the above steps 8 & 9 are repeated for tmle_g2_out and the ATE.
    # First restore (Anodes,sA) that were generated in the observed bootstrapped sample!
    # In this case the intervention summaries, s.a, "def_new_sA(A = A)" will correctly use the observed bootstrapped values of A
    # -----------------------------------------------------------------------------------------------------------------------------
    if (!is.null(tmle_g2_out)) {
      # * restore the observed boostrapped Anodes and sA
      DatNet.ObsP0$Odata$restoreAnodes()
      # * verify sA's were also restored and if not, regenerate them
      if (!DatNet.ObsP0$Odata$restored_sA_Vars) DatNet.ObsP0$datnetA$make.sVar(sVar.object = tmle_g2_out$sA)
      DatNet.gstar$make.dat.sWsA(p = 1, new.sA.object = tmle_g2_out$new.sA, sA.object = tmle_g2_out$sA, DatNet.ObsP0 = DatNet.ObsP0)

      fWi.g2 <- psi.evaluator$get.gcomp(m.Q.init = tmle_g2_out$m.Q.init)
      tmle.g2 <- psi.evaluator$get.tmleB(m.Q.starB = tmle_g2_out$MC_fit_params$m.Q.star)
      mean_obs_gcomp_g2 <- mean_obs_gcomp_g2 + fWi.g2
      mean_obs_tmle_B_g2 <- mean_obs_tmle_B_g2 + tmle.g2
    }
  }

  mean_obs_gcomp_g1 <- mean_obs_gcomp_g1 / n.MC
  mean_obs_gcomp_g2 <- mean_obs_gcomp_g2 / n.MC
  mean_obs_tmle_B_g1 <- mean_obs_tmle_B_g1 / n.MC
  mean_obs_tmle_B_g2 <- mean_obs_tmle_B_g2 / n.MC

  # Restore the original data.table OdataDT, Y values and indicator of deterministic Y values
  # Otherwise it messes up the IC-based inference, since it uses the observed Yvals to est. the IC
  DatNet.ObsP0$Odata$OdataDT <- OdataDT.P0
  DatNet.ObsP0$noNA.Ynodevals <- noNA.Ynodevals.P0
  DatNet.ObsP0$det.Y <- det.Y.P0

  # mean_obs_gcomp_g1 = mean_obs_gcomp_g1,
  out_mean_tmleB <- list(EY_gstar1 = mean_obs_tmle_B_g1, EY_gstar2 = mean_obs_tmle_B_g2, ATE = mean_obs_tmle_B_g1 - mean_obs_tmle_B_g2)
  return(out_mean_tmleB)
}



# --------------------------------------------------------------------------------------------------------
# Bootstrap with resampling of (sW,sA,Y) with replacement (as if iid) - DOESN'T WORK FOR DEPENDENT DATA
# --------------------------------------------------------------------------------------------------------
iid_bootstrap_tmle <- function(n.boot, estnames, DatNet.ObsP0, tmle_g_out, QY_mat, wts_mat) {
  QY.init <- QY_mat[, "QY.init"]
  off <- qlogis(QY.init)  # offset
  DatNet.gstar <- tmle_g_out$DatNet.gstar
  psi.evaluator <- tmle_g_out$psi.evaluator
  m.Q.init <- tmle_g_out$m.Q.init
  m.h.fit <- tmle_g_out$m.h.fit
  h_wts <- wts_mat[, "h_wts"]
  Y <- DatNet.ObsP0$noNA.Ynodevals
  determ.Q <- DatNet.ObsP0$det.Y

  #************************************************
  # BOOTSTRAP TMLE UPDATES:
  #************************************************
  # n.boot <- 100
  boot_eps <- vector(mode = "numeric", length = n.boot)
  boot_tmle_B <- vector(mode = "numeric", length = n.boot)
  for (i in (1:n.boot)) {
      boot_idx <- sample.int(n = DatNet.ObsP0$nobs, replace = TRUE)
      boot.tmle.obj <- tmle.update(estnames = estnames,
                                   Y = Y[boot_idx], off = off[boot_idx], h_wts = h_wts[boot_idx],
                                   determ.Q = determ.Q[boot_idx], predictQ = FALSE)
      boot_eps[i] <- boot.tmle.obj$m.Q.star.coef
      boot_tmle_B[i] <- mean(psi.evaluator$get.boot.tmleB(m.Q.starB = boot_eps[i], boot_idx = boot_idx))
  }
  var_tmleB_boot <- var(boot_tmle_B)

  # print("mean(boot_eps): "); print(mean(boot_eps))
  # print("mean(boot_tmle_B): "); print(mean(boot_tmle_B))
  # print("var(boot_tmle_B): "); print(var(boot_tmle_B))
  return(var_tmleB_boot)
}

# --------------------------------------------------------------------------------------------------------
# Parametric bootstrap: sampling W as iid, boot.nodes from m.h.fit or special fit obj and Y from m.Q.init.N;
# Parametric bootstrap: need to explicitly spec. fitted nodes which are cond independent and then resample from those fits.
# --------------------------------------------------------------------------------------------------------
# Could be very different from Anodes or part of Anodes. If no reg form was specified, finds it among hforms by default;
# However, Y's are always re-sampled from m.Q.init
# Allows boot.nodes to be different from Anodes with their own regression formulas
# These regressions are specified with "boot.form" arg: specifies regression forms for conditionally independent nodes to be resampled from such fits
par_bootstrap_tmle <- function(n.boot, boot.nodes, boot.form, estnames, DatNet.ObsP0, tmle_g1_out, tmle_g2_out) {
  # ******** REPLACED this with the actual model fit g.N *********
  f.g0 <- NULL

  boot_eps_g1 <- vector(mode = "numeric", length = n.boot)
  boot_eps_g2 <- vector(mode = "numeric", length = n.boot)

  boot_gcomp_g1 <- vector(mode = "numeric", length = n.boot)
  boot_gcomp_g2 <- vector(mode = "numeric", length = n.boot)
  boot_gcomp_ATE <- vector(mode = "numeric", length = n.boot)

  boot_tmle_B_g1 <- vector(mode = "numeric", length = n.boot)
  boot_tmle_B_g2 <- vector(mode = "numeric", length = n.boot)
  boot_tmle_B_ATE <- vector(mode = "numeric", length = n.boot)

  boot_IC_tmle <- vector(mode = "numeric", length = n.boot)

  psi.evaluator <- tmle_g1_out$psi.evaluator
  # Always start with the observed Anodes
  if (!DatNet.ObsP0$Odata$curr_data_A_g0) {
    DatNet.ObsP0$Odata$restoreAnodes()
  }

  # -----------------------------------------------------------------------------------------------
  # Save the original input data.table OdataDT, otherwise it will be over-written:
  # Note we are not saving the original saved Anodes and sA -> these fields NULLed during bootstrap
  # -----------------------------------------------------------------------------------------------
  OdataDT.P0 <- DatNet.ObsP0$Odata$OdataDT
  noNA.Ynodevals.P0 <- DatNet.ObsP0$noNA.Ynodevals
  det.Y.P0 <- DatNet.ObsP0$det.Y
  DatNet.gstar <- tmle_g1_out$DatNet.gstar

  # -----------------------------------------------------------------------------------------------
  # fit any regressions specified in boot.form:
  # -----------------------------------------------------------------------------------------------
  if (!is.null(boot.form)) {
    boot.form_fit <- fit_reg.forms(boot.form, DatNet.ObsP0, tmle_g1_out$sW, tmle_g1_out$sA)
    summeas_boot <- boot.form_fit$summeas_boot
    boot_outvar_nms <- unlist(boot.form_fit$boot_outvar_nms)
  } else {
    boot_outvar_nms <- NULL
  }

  # -----------------------------------------------------------------------------------------------
  # loop over n.boot
  # -----------------------------------------------------------------------------------------------
  for (i in (1:n.boot)) {
    # 1. Resample W (with replacement) by re-purposing the instance of DatNet.gstar; Re-evaluate baseline summaries sW; Re-shuffle pre-saved values of Y and det.Y;
    boot_idx <- sample.int(n = DatNet.ObsP0$nobs, replace = TRUE)
    DatNet.ObsP0$Odata$OdataDT <- OdataDT.P0[boot_idx, ]

    # Need to NULL previously backed-up values of A and sA, since they no longer correspond with bootstrapped sample
    DatNet.ObsP0$Odata$A_g0_DT <- NULL
    DatNet.ObsP0$Odata$sA_g0_DT <- NULL
    DatNet.ObsP0$Odata$save_sA_Vars <- NULL

    DatNet.ObsP0$datnetW$make.sVar(sVar.object = tmle_g1_out$sW)
    DatNet.ObsP0$datnetW$fixmiss_sVar() # permanently replace NA values in sW with 0
    DatNet.ObsP0$det.Y <- det.Y.P0[boot_idx]
    DatNet.ObsP0$noNA.Ynodevals <- noNA.Ynodevals.P0[boot_idx]

    # 2. Generate new A's from g0 or g.N (replace A with sampled A's in DatNet.ObsP0); Re-evaluate exposure summaries (sA) based on new DatNet.ObsP0:
    if (is.null(f.g0)) {
      for (Anode in boot.nodes) {
        # If Anode was referenced in boot.form, then sample Anode from this new fit, otherwise sample from hforms fit
        if (Anode %in% boot_outvar_nms) {
          model.sVar.gN <- findRegSummaryObj(summeas_boot, outvar = Anode)
        } else {
          # NOTE: new RegressionClass implementation will allow pulling top level outvar names directly
          # (regardless of how the regressions were specified ("A1+A2~W" or c("A1~W","A2~W")))
          # will allow checking directly if model for Anode exists: (Anode %in% unlist(tmle_g1_out$m.h.fit$summeas.g0$reg$outvar))
          # current implementation doesn't allow that, so need to check the model is not NULL
          model.sVar.gN <- findRegSummaryObj(tmle_g1_out$m.h.fit$summeas.g0, outvar = Anode)
        }
        if (is.list(model.sVar.gN)) model.sVar.gN <- model.sVar.gN[[1]]
        if (is.null(model.sVar.gN)) stop("Model fit for parametric bootstrap variable '" %+% Anode %+% "' could not be located. " %+%
          "Make sure to specify the regression formula for every variable in 'boot.nodes', either as part of regression formulas in 'hform.g0'/'hform.gstar'" %+%
          "or 'boot.form'")
        A.sample.gN <- model.sVar.gN$sampleA(newdata = DatNet.ObsP0)
        DatNet.ObsP0$Odata$replaceOneAnode(AnodeName = Anode, newAnodeVal = A.sample.gN)
      }
    } else if (!is.null(f.g0)) {
      DatNet.ObsP0$make.dat.sWsA(p = 1, f.g_fun = f.g0, sA.object = tmle_g1_out$sA, DatNet.ObsP0 = DatNet.ObsP0)
    }

    DatNet.ObsP0$datnetA$make.sVar(sVar.object = tmle_g1_out$sA) # re-evaluate exposure summaries based on bootstraped W and re-sampled A's
    DatNet.ObsP0$Odata$curr_data_A_g0 <- TRUE

    # 4. Predict P(Y_i=1|sW,sA) using m.Q.init (the initial fit \bar{Q}_N) based on newly resampled (sW,sA):
    detY.boot <- DatNet.ObsP0$det.Y
    QY.init.boot <- DatNet.ObsP0$noNA.Ynodevals
    QY.init.boot[!detY.boot] <- tmle_g1_out$m.Q.init$predict(newdata = DatNet.ObsP0)$getprobA1[!detY.boot] # getting predictions P(Y=1) for non-DET Y
    off.boot <- qlogis(QY.init.boot)  # offset

    # 5. Sample a vector of new (Y_i, i=1,...,N):
    # **** TO BE REPLACED WITH SOMETHING LIKE B psi.evaluator$sampleY(Qprob = QY.init) ****
    Y.boot <- rbinom(n = length(QY.init.boot), size = 1, prob = QY.init.boot)

    # 6. Predict new weights h_wts = P_{\bar{g}^*_N}(sA | sW)/P_{\bar{g}_0}(sA | sW) :
    # using previously fitted m.h.fit$summeas.g0 and m.h.fit$summeas.gstar and the newly resampled (sW,sA) under Q.W.N, m.g.N
    h_wts_g1.boot <- predict.hbars(newdatnet = DatNet.ObsP0, m.h.fit = tmle_g1_out$m.h.fit)

    # 7. Fit a TMLE update epsilon on this new bootstapped dataset.
    boot.tmle_g1.obj <- tmle.update(estnames = estnames,
                                    Y = Y.boot, off = off.boot, h_wts = h_wts_g1.boot,
                                    determ.Q = detY.boot, predictQ = FALSE)
    boot_eps_g1[i] <- boot.tmle_g1.obj$m.Q.star.coef

    # If tmle_g2_out is present, do the same, before DatNet.ObsP0 is overwritten
    if (!is.null(tmle_g2_out)) {
      h_wts_g2.boot <- predict.hbars(newdatnet = DatNet.ObsP0, m.h.fit = tmle_g2_out$m.h.fit)
      boot.tmle_g2.obj <- tmle.update(estnames = estnames,
                                      Y = Y.boot, off = off.boot, h_wts = h_wts_g2.boot,
                                      determ.Q = detY.boot, predictQ = FALSE)
      boot_eps_g2[i] <- boot.tmle_g2.obj$m.Q.star.coef
    }

    # 8. Re-create DatNet.gstar with boostrapped summaires sW and sA generated under f.gstar:
    # However, first need to save (Anodes,sA) from the observed bootstrapped sample, otherwise they are forever lost
    DatNet.ObsP0$Odata$backupAnodes(sA = tmle_g1_out$sA)
    # Will over-write Anodes/sA in DatNet.ObsP0$Odata:
    DatNet.gstar$make.dat.sWsA(p = 1, new.sA.object = tmle_g1_out$new.sA, sA.object = tmle_g1_out$sA, DatNet.ObsP0 = DatNet.ObsP0)

    # 9. Evaluate the substitution estimator and the components of the EIC D_Y and D_W:
    fWi.boot_g1 <- psi.evaluator$get.gcomp(m.Q.init = tmle_g1_out$m.Q.init)
    boot_gcomp_g1[i] <- mean(fWi.boot_g1)
    boot_tmle_B_g1[i] <- mean(psi.evaluator$get.tmleB(m.Q.starB = boot.tmle_g1.obj$m.Q.star.coef))

    obsYvals.boot = DatNet.ObsP0$noNA.Ynodevals
    boot_IC_tmle[i] <- mean(h_wts_g1.boot * (obsYvals.boot - QY.init.boot) + (fWi.boot_g1 - boot_tmle_B_g1[i]))

    # 10. If (!is.null(tmle_g2_out)) then the above steps 8 & 9 are repeated for tmle_g2_out and the ATE.
    # First restore (Anodes,sA) that were generated in the observed bootstrapped sample!
    # In this case the intervention summaries, s.a, "def_new_sA(A = A)" will correctly use the observed bootstrapped values of A
    # -----------------------------------------------------------------------------------------------------------------------------
    # Q: How do we get the bootsrap var for the ATE?
      # 1) We can generate a new TMLE update that directly corresponds to the EIC for ATE (need to figure out EIC)
      # 2) We could also evaluate the ATE as a plug-in estimator from two separate TMLE updates - using this approach.
    # Both approaches should be equivalent, except for situations with 2 "incompatible" target parameters
    # -----------------------------------------------------------------------------------------------------------------------------
    if (!is.null(tmle_g2_out)) {
      # * restore the observed boostrapped Anodes and sA
      DatNet.ObsP0$Odata$restoreAnodes()
      # * verify sA's were also restored and if not, regenerate them
      if (!DatNet.ObsP0$Odata$restored_sA_Vars) DatNet.ObsP0$datnetA$make.sVar(sVar.object = tmle_g2_out$sA)
      DatNet.gstar$make.dat.sWsA(p = 1, new.sA.object = tmle_g2_out$new.sA, sA.object = tmle_g2_out$sA, DatNet.ObsP0 = DatNet.ObsP0)

      fWi.boot_g2 <- psi.evaluator$get.gcomp(m.Q.init = tmle_g2_out$m.Q.init)
      boot_gcomp_g2[i] <- mean(fWi.boot_g2)
      boot_gcomp_ATE[i] <- boot_gcomp_g1[i] - boot_gcomp_g2[i]

      boot_tmle_B_g2[i] <- mean(psi.evaluator$get.tmleB(m.Q.starB = boot.tmle_g2.obj$m.Q.star.coef))
      boot_tmle_B_ATE[i] <- boot_tmle_B_g1[i] - boot_tmle_B_g2[i]
    }
  }

  var_gcomp_boot_g1 <- var(boot_gcomp_g1)
  var_gcomp_boot_g2 <- var(boot_gcomp_g2)
  var_gcomp_boot_ATE <- var(boot_gcomp_ATE)

  var_tmleB_boot_g1 <- var(boot_tmle_B_g1)
  var_tmleB_boot_g2 <- var(boot_tmle_B_g2)
  var_tmleB_boot_ATE <- var(boot_tmle_B_ATE)

  # browser()
  # boxplot(nF.PA.tab_mat[,c(1:5)], main = "bootstrap distr of nF.PA over 500 reps")
  # boxplot(nF.PA.tab_mat[,c(1:10)], main = "bootstrap distr of nF.PA over 500 reps")
  # plot(density(boot_gcomp_g1))
  # plot(density(boot_tmle_B_g1))

  # print("boot_time for n.boot = " %+% n.boot); print(boot_time)
  # [1] "boot_time for n.boot = 500"
  #    user  system elapsed
  # 708.094 134.648 844.491

  if (!is.null(tmle_g2_out)) {
    gcomp_g1_boot_col <- rbind(tmle_g1_out$ests_mat["MLE",],                                  mean(boot_gcomp_g1), var(boot_gcomp_g1))
    gcomp_g2_boot_col <- rbind(tmle_g2_out$ests_mat["MLE",],                                  mean(boot_gcomp_g2), var(boot_gcomp_g2))
    gcomp_ATE_boot_col <- rbind(tmle_g1_out$ests_mat["MLE",]-tmle_g2_out$ests_mat["MLE",],    mean(boot_gcomp_ATE), var(boot_gcomp_ATE))
    tmle_B_g1_boot_col <- rbind(tmle_g1_out$ests_mat["TMLE",],                                mean(boot_tmle_B_g1), var(boot_tmle_B_g1))
    tmle_B_g2_boot_col <- rbind(tmle_g2_out$ests_mat["TMLE",],                                mean(boot_tmle_B_g2), var(boot_tmle_B_g2))
    tmle_B_ATE_boot_col <- rbind(tmle_g1_out$ests_mat["TMLE",]-tmle_g2_out$ests_mat["TMLE",], mean(boot_tmle_B_ATE), var(boot_tmle_B_ATE))

    res_mat <- cbind(gcomp_g1 = gcomp_g1_boot_col, gcomp_g2 = gcomp_g2_boot_col, gcomp_ATE = gcomp_ATE_boot_col,
                     tmle_B_g1 = tmle_B_g1_boot_col, tmle_B_g2 = tmle_B_g2_boot_col, tmle_B_ATE = tmle_B_ATE_boot_col)

    rownames(res_mat) <- c("Mean.Est", "Boot.Mean.Est", "Boot.Var")
    colnames(res_mat) <- c("gcomp_g1","gcomp_g2","gcomp_ATE","tmle_B_g1","tmle_B_g2","tmle_B_ATE")
    print("Parametric Bootstrap Inference: "); print(res_mat)

    est_mat <- cbind(tmle_g1_out$ests_mat, tmle_g2_out$ests_mat, tmle_g1_out$ests_mat - tmle_g2_out$ests_mat)
    colnames(est_mat) <- c("g1", "g2", "ATE")
    print("est_mat"); print(est_mat)
  }

  # Restore the original data.table OdataDT, Y values and indicator of deterministic Y values
  # Otherwise it messes up the IC-based inference, since it uses the observed Yvals to est. the IC
  DatNet.ObsP0$Odata$OdataDT <- OdataDT.P0
  DatNet.ObsP0$noNA.Ynodevals <- noNA.Ynodevals.P0
  DatNet.ObsP0$det.Y <- det.Y.P0

  out_var_tmleB_boot <- list(EY_gstar1 = var_tmleB_boot_g1, EY_gstar2 = var_tmleB_boot_g2, ATE = var_tmleB_boot_ATE)
  return(out_var_tmleB_boot)
}

# create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
# boot.var, n.boot,
make_EYg_obj <- function(estnames, estoutnames, alpha, DatNet.ObsP0, tmle_g_out, tmle_g2_out=NULL, MC.tmle.eval = NULL, var_tmleB_boot) {
  nobs <- DatNet.ObsP0$nobs
  NetInd_k <- DatNet.ObsP0$netind_cl$NetInd_k
  nF <- DatNet.ObsP0$netind_cl$nF
  ests_mat <- tmle_g_out$ests_mat
  QY_mat <- tmle_g_out$QY_mat
  fWi_mat <- tmle_g_out$fWi_mat
  wts_mat <- tmle_g_out$wts_mat

  boot.as.var_tmleB <- var_tmleB_boot * nobs # parametric bootstrap-based asymptotic variance (est. outside this function)

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
                              MC.tmle.eval = MC.tmle.eval,
                              ests_mat = ests_mat,
                              QY_mat = QY_mat,
                              wts_mat = wts_mat, fWi_mat = fWi_mat)
  )
  if (gvars$verbose) {
    print("time to estimate Vars: "); print(getVar_time)
  }

  # ------------------------------------------------------------------------------------------
  # parametric bootstrap-based asymptotic variance estimates matrix:
  # ------------------------------------------------------------------------------------------
  boot.as.var_mat <- matrix(nrow = nrow(as.vars_obj$as.var_mat), ncol = 1)
  boot.as.var_mat[1,1] <- boot.as.var_tmleB
  rownames(boot.as.var_mat) <- rownames(as.vars_obj$as.var_mat)
  colnames(boot.as.var_mat) <- colnames(as.vars_obj$as.var_mat)

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

  # CIs based on bootstrapped variance:
  boot.CIs_mat <- t(apply(cbind(ests_mat, boot.as.var_mat), 1, get_CI, n = nobs))
  colnames(boot.CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))

  # CIs based on the IC for dependent data:
  CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$as.var_mat), 1, get_CI, n = nobs))
  colnames(CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))

  # CIs based on conditional on W variance:
  condW.CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$condW.vars_mat), 1, get_CI, n = nobs))
  colnames(condW.CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))

  # CIs based on IID variance:
  iid.CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$iid.vars_mat), 1, get_CI, n = nobs))
  colnames(iid.CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))

  # ------------------------------------------------------------------------------------------
  # Rename estimators for the final output:
  # ------------------------------------------------------------------------------------------
  rownames(ests_mat) <- estoutnames
  rownames(as.vars_obj$as.var_mat) <- rownames(boot.as.var_mat) <- rownames(as.vars_obj$condW.vars_mat) <- rownames(as.vars_obj$iid.vars_mat) <- estoutnames
  rownames(CIs_mat) <- rownames(boot.CIs_mat) <- rownames(condW.CIs_mat) <- rownames(iid.CIs_mat) <- estoutnames

  EY_g.star <- list(estimates = ests_mat,

                    boot.vars = (boot.as.var_mat / nobs), # parametric bootstrap variance
                    IC.vars = (as.vars_obj$as.var_mat / nobs), # IC-based variance (dependent)
                    condW.IC.vars = (as.vars_obj$condW.vars_mat / nobs), # IC-based variance conditional on W
                    iid.vars = (as.vars_obj$iid.vars_mat / nobs), # iid Variance

                    aux.vars = (as.vars_obj$aux.vars_mat / nobs), # auxilary (additional variance estimators)

                    boot.CIs = boot.CIs_mat,
                    IC.CIs = CIs_mat,
                    condW.CIs = condW.CIs_mat,
                    iid.CIs = iid.CIs_mat,

                    h_g0_SummariesModel = NULL,
                    h_gstar_SummariesModel = NULL)

  if (is.null(tmle_g2_out)) {
    EY_g.star[["h_g0_SummariesModel"]] <- tmle_g_out$h_g0_SummariesModel
    EY_g.star[["h_gstar_SummariesModel"]] <- tmle_g_out$h_gstar_SummariesModel
  }

  return(EY_g.star)
}