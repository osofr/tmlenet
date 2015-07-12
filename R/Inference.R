  #----------------------------------------------------------------------------------
  # Estimate IC-based Variance (as. Var) and CIs (based on f_W and fY)
  #----------------------------------------------------------------------------------
  # use estimates of fWi (hold all W's fixed at once),
  # loop over all intersecting friends networks
  # calculate R_ij*(fW_i*fW_j) - see page 33 Network paper vdL
  #----------------------------------------------------------------------------------
  # Helper function to calculate cross product sum of correlated f_Wi (see p.33 of vdL)
  # New fast method for as var calculation (matrix vs)
  est_sigmas <- function(n, NetInd_k, obsYvals, ests_mat, QY_mat, wts_mat, fWi_mat, onlyTMLE_B) {
      fWi <- fWi_mat[, "fWi_Qinit"]
      QY.init <- QY_mat[, "QY.init"] # QY.star <- QY_mat[, "QY.star_A"]
      h_wts <- wts_mat[,"h_wts"]
      g_wts <- wts_mat[,"g_wts"]

      var_tmle_A <- var_tmleiptw_1stO <- var_tmleiptw_2ndO <- 0
      var_tmle_A_Q.init <- var_tmle_B_Q.init <- 0

      # get the connectivity n_by_n mtx (1 indicates intersection of friendship sets)
      # returns 1) 1st order and 2) 1st and 2nd order connections
      get.Fiintersectmtx <- function(n) {
          # turn NetInd_k into an n_by_n mtx of indicators (friends=1 vs not friends=0)
          # this is not symmetric, friends are rows, connections (friends set) are cols
          # obs i is always in i's friend set Fi
          # expand the friend network by including friends of friends
          # to get the second-order connections (intersections of friends of friends)
          NetInd_k_2ndO <- NetInd_k
          for (j in (1:ncol(NetInd_k))) {
            NetInd_k_2ndO <- cbind(NetInd_k_2ndO, NetInd_k[NetInd_k[,j],])
          }
          conn_ind_mtx_1st <- diag(x=1,n,n)
          conn_ind_mtx_2ndO <- diag(x=1,n,n)
          for (i in (1:n)) { # CAN REPLACE WITH A LOOP OVER COLUMNS of NetInd_k INSTEAD!!!!!!!!
              conn_ind_mtx_1st[i, NetInd_k[i,]] <- 1 # non-symmetric connectivity mtx
              conn_ind_mtx_2ndO[i, NetInd_k_2ndO[i,]] <- 1 # capture second order connections as well (for g_i and g_j that are dependent)
          }
          # second pass over columns of connectivity mtx to connect indirect intersections (i and j have a common friend but are not friends themselves):
          conn_ind_mtx_1st_indir <- conn_ind_mtx_1st
          conn_ind_mtx_2ndO_indir <- conn_ind_mtx_2ndO
          for (j in (1:n)) {
            conn_ind_mtx_1st_indir[which(conn_ind_mtx_1st[, j]==1), which(conn_ind_mtx_1st[, j]==1)] <- 1
            conn_ind_mtx_2ndO_indir[which(conn_ind_mtx_2ndO[, j]==1), which(conn_ind_mtx_2ndO[, j]==1)] <- 1
          }
          return(list(conn_ind_mtx_1st=conn_ind_mtx_1st_indir, conn_ind_mtx_2ndO=conn_ind_mtx_2ndO_indir))
      }
      get.crossprodmtx <- function(fvec_i) { # cross prod mtx (n_by_n) for any vector of size n
          return(tcrossprod(fvec_i, fvec_i))
      }
      est.sigma_fsum <- function(fcrossprod, connectmtx)  { # sum the cross prod vector over connectivity mtx (prod will only appear if (i,j) entry is 1)
          return((1/n) * sum(fcrossprod * connectmtx))
      }
      ###################################
      # Capturing dependencies for IPTW-based estimators (second-order dependencies)
      # Algorithm:
      # Extend the friends network for each i to reflect that for i, g(A|W) is a product over j\inF_i of g(A_j|W_m:m\inF_j)
      # Thus, create a new set of F_star_i for i, where F_star_i = union(F_j:j\inF_i)
      ###################################
      connectmtx_obj <- get.Fiintersectmtx(n = n)
      connectmtx_1stO <- connectmtx_obj$conn_ind_mtx_1st
      connectmtx_2ndO <- connectmtx_obj$conn_ind_mtx_2ndO

      # TMLE A (clever covariate update): Inference based on the iid IC analogy, QY.init := initial Q model predictions, h_wts := h_tilde
      if (!onlyTMLE_B) {
        # iidIC_tmle_A <- h_wts * (obsYvals - QY.init) + fWi_A
        iidIC_tmle_A <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",])
        var_tmle_A <- est.sigma_fsum(get.crossprodmtx(iidIC_tmle_A), connectmtx_1stO)
      }

      # TMLE B (weighted model update): Inference based on the iid IC
      # iidIC_tmle_B <- h_wts * (obsYvals - QY.init) + fWi_B
      iidIC_tmle_B <- h_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",])
      var_tmle_B <- est.sigma_fsum(get.crossprodmtx(iidIC_tmle_B), connectmtx_1stO)
      # simple iid estimator of the asymptotic variance (no adjustment made when two observations i!=j are dependent):
      var_iid.tmle_B <- mean((iidIC_tmle_B)^2)

      # TMLE based on iptw clever covariate (more non-parametric)
      if (!onlyTMLE_B) {
        iidIC_tmleiptw <- g_wts * (obsYvals - QY.init) + (fWi - ests_mat[rownames(ests_mat)%in%"tmle_g_iptw",])
        var_tmleiptw_1stO <- est.sigma_fsum(get.crossprodmtx(iidIC_tmleiptw), connectmtx_1stO)
        var_tmleiptw_2ndO <- est.sigma_fsum(get.crossprodmtx(iidIC_tmleiptw), connectmtx_2ndO)
      }

      # IPTW h (based on the mixture density clever covariate (h)):
      iidIC_iptw_h <- h_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"h_iptw",])
      var_iptw_h <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw_h), connectmtx_1stO)

      # IPTW g:
      iidIC_iptw_g <- g_wts * (obsYvals) - (ests_mat[rownames(ests_mat)%in%"g_iptw",])
      var_iptw_1stO <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw_g), connectmtx_1stO)
      var_iptw_2ndO <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw_g), connectmtx_2ndO)
 
      # Inference based on the EIC, with factorization into orthogonal components sigma2_DY and sigma2_W_N
      # sigma2_DY_i are independent (since they are conditioned on W,A)
      # sigma2_W_N_i are dependent => need to take double sum of their crossprod among dependent units
      if (!onlyTMLE_B) {
        D_star_Yi.Qinit <- h_wts * (obsYvals - QY.init) # h*(Y-Q_bar_N):
        sigma2_DY <- (1/n) * sum(D_star_Yi.Qinit^2)  # Sum_{i} (D_star_Yi)^2

        fW_A_crossprod <- get.crossprodmtx((fWi - ests_mat[rownames(ests_mat)%in%"tmle_A",]))
        sigma2_W_N_A <- est.sigma_fsum(fW_A_crossprod, connectmtx_1stO)
        var_tmle_A_Q.init <- sigma2_W_N_A + sigma2_DY

        # **NEW** TMLE B (weights model update)
        fW_B_crossprod <- get.crossprodmtx((fWi - ests_mat[rownames(ests_mat)%in%"tmle_B",]))
        sigma2_W_N_B <- est.sigma_fsum(fW_B_crossprod, connectmtx_1stO)
        var_tmle_B_Q.init <- sigma2_W_N_B + sigma2_DY

        # D_star_Yi.Qstar <- h_wts * (obsYvals - QY.star)
        # D_star_Yi.Qstar[determ.Q] <- 0
        # fDY_crossprod <- get.crossprodmtx(D_star_Yi.Qstar)
        # double sum over dependent subjects, Sum_{i,j} R_W(i,j)*D_star_Yi*D_star_Yj
        # sigma2_Y_N <- est.sigma_fsum(fDY_crossprod, connectmtx_1stO)
        # var_tmle_Q.init_c <- sigma2_W_N_A + sigma2_Y_N

        # #--------
        # # conservative estimate of the as. variance from EIC for TMLE A:
        # # abs terms double sum over dependent subjects, Sum_{i,j} R_W(i,j)*|D_star_Yi|*|D_star_Yj|:
        # fabsDY_crossprod <- get.crossprodmtx(abs(D_star_Yi.Qstar))
        # abs_sigma2_Y_N <- est.sigma_fsum(fabsDY_crossprod, connectmtx_1stO)
        # var_tmle_A_Q.star_cons <- sigma2_W_N_A + abs_sigma2_Y_N
        # # --------
      }

      var.ests <- c(abs(var_tmle_A), abs(var_tmle_B), abs(var_tmleiptw_1stO), abs(var_iptw_h), abs(var_iptw_1stO), 0)
      estnames <- c(      "tmle_A",     "tmle_B",       "tmle_g_iptw",        "h_iptw",         "g_iptw", "mle")
      as.var_mat <- matrix(0, nrow = length(var.ests), ncol = 1)
      as.var_mat[,1] <- var.ests
      rownames(as.var_mat) <- estnames
      colnames(as.var_mat) <- "var"

      other.vars = c(
                    var_iid.tmle_B = abs(var_iid.tmle_B), # no adjustment for correlations i,j
                    var_tmleiptw_2ndO = abs(var_tmleiptw_2ndO), # adjusting for 2nd order dependence of i,j
                    var_iptw_2ndO = abs(var_iptw_2ndO), # adjusting for 2nd order dependence of i,j
                    var_tmle_A_Q.init = abs(var_tmle_A_Q.init), # using the EIC & Q.init for TMLE A
                    var_tmle_B_Q.init = abs(var_tmle_B_Q.init)  # using the EIC & Q.init for TMLE B
                    )

      return(list(as.var_mat = as.var_mat, other.vars = other.vars))
                  # other.vars = list(
                  #   var_iid.tmle_B = abs(var_iid.tmle_B), # no adjustment for correlations i,j
                  #   var_tmleiptw_2ndO = abs(var_tmleiptw_2ndO), # adjusting for 2nd order dependence of i,j
                  #   var_iptw_2ndO = abs(var_iptw_2ndO), # adjusting for 2nd order dependence of i,j
                  #   var_tmle_A_Q.init = abs(var_tmle_A_Q.init), # using the EIC & Q.init for TMLE A
                  #   var_tmle_B_Q.init = abs(var_tmle_B_Q.init)  # using the EIC & Q.init for TMLE B
                  # )
          # ))
  }

# create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
make_EYg_obj <- function(alpha, onlyTMLE_B, datNetObs, tmle_g_out, tmle_g2_out=NULL) {
  # print("tmle_g_out$ests_mat[, tmle_A]"); print(tmle_g_out$ests_mat[rownames(tmle_g_out$ests_mat)%in%"tmle_A",])
  nobs <- datNetObs$nobs
  NetInd_k <- datNetObs$netind_cl$NetInd_k

  ests_mat <- tmle_g_out$ests_mat
  QY_mat <- tmle_g_out$QY_mat
  fWi_mat <- tmle_g_out$fWi_mat
  wts_mat <- tmle_g_out$wts_mat
  if (!is.null(tmle_g2_out)) {
    ests_mat <- tmle_g_out$ests_mat - tmle_g2_out$ests_mat
    fWi_mat <- tmle_g_out$fWi_mat - tmle_g2_out$fWi_mat
    wts_mat <- tmle_g_out$wts_mat - tmle_g2_out$wts_mat
  }
  # get the iid IC-based asymptotic variance estimates:
  as.vars_obj <- est_sigmas(n = nobs, NetInd_k = NetInd_k, obsYvals = datNetObs$noNA.Ynodevals,
                            ests_mat = ests_mat, QY_mat = QY_mat, wts_mat = wts_mat, fWi_mat = fWi_mat,
                            onlyTMLE_B = onlyTMLE_B)

 get_CI <- function(xrow, n) {
    f_est_CI <- function(n, psi, sigma2_N) { # get CI
      z_alpha <- qnorm(1-alpha/2)
      CI_est <- c(psi - z_alpha*sqrt(sigma2_N) / sqrt(n), psi + z_alpha*sqrt(sigma2_N) / sqrt(n))
      return(CI_est)
    }
    psi <- xrow["estimate"];
    sigma2_N <- xrow["var"];
    return(f_est_CI(n = n, psi = psi, sigma2_N = sigma2_N))
  }

  CIs_mat <- t(apply(cbind(ests_mat, as.vars_obj$as.var_mat), 1, get_CI, n = nobs))
  colnames(CIs_mat) <- c("LBCI_"%+%as.character(alpha/2), "UBCI_"%+%as.character(1-alpha/2))

  EY_g.star <- list(estimates = ests_mat,
                    vars = (as.vars_obj$as.var_mat / nobs),
                    CIs = CIs_mat,
                    other.vars = (as.vars_obj$other.vars / nobs)
                    # other.vars = lapply(as.vars_obj$other.vars, function(var) var / nobs)
                    )

  return(EY_g.star)
}