
# (DEPRECATED)
#-----------------------------------------------------------------------------
# Fit glm based on formula in "form"
#-----------------------------------------------------------------------------
f_est <- function(d, form, family) {
  ctrl <- glm.control(trace=FALSE, maxit=1000)
    SuppressGivenWarnings({
              m <- glm(as.formula(form),
                  data=d,
                  family=family,
                  control=ctrl)
              },
              GetWarningsToSuppress())
    return(m)
}


# (DEPRECATED)
#-----------------------------------------------------------------------------
# USE glm.fit FUNCTION FOR FASTER FITTING of LOGISTIC REG TAKES DESIGN MAT AND Y VECTOR
#-----------------------------------------------------------------------------
.f.est_binom_fast <- function(X_mat, Y_vals) {
  ctrl <- glm.control(trace=FALSE, maxit=1000)          
    SuppressGivenWarnings({
              m.fit <- glm.fit(x=cbind(1,X_mat), y=Y_vals, family = binomial(), control=ctrl)
              }, GetWarningsToSuppress())
    return(m.fit)
}
.f_predict_fast <- function(glmfit, new_mtx) {
      new_mtx <- cbind(1,new_mtx)
      eta <- new_mtx[,!is.na(glmfit$coef), drop=FALSE] %*% glmfit$coef[!is.na(glmfit$coef)]
      return(glmfit$family$linkinv(eta))
}

# (DEPRECATED) (WAS USED FOR PLUG-IN NPMLE ESTIMATOR OF \bar{h}_gN)
# sample treatments with probA from above fcn
.f.gen.A_N <- function(df, deterministic, m.gN) {
	n <- nrow(df)
  	rbinom(n, 1, .f.gen.probA_N(df, deterministic, m.gN))
}


# (DEPRECATED)
#-----------------------------------------------------------------------------
# Calculate the joint probability of observing a=(a_1,..,a_k) from P(A^s=1)
# takes the matrix of predicted probabilities P(A^s=1): data.probA
# and a matrix of observed a^s (a_1,..,a_k): data.indA
#-----------------------------------------------------------------------------
.f.cumprod.matrix <- function(data.indA, data.probA) {
  y <- matrix(1, nrow=dim(data.probA)[1], ncol=dim(data.probA)[2])
  y[, 1] <- data.probA[,1]^as.integer(data.indA[,1]) *
              (1-data.probA[,1])^(1-as.integer(data.indA[,1]))
  if (dim(data.probA)[2] > 1) {
    for (i in 2:dim(data.probA)[2]) {
      y[,i] <- y[,i-1] * (data.probA[,i]^as.integer(data.indA[,i]) * 
                (1-data.probA[,i])^(1-as.integer(data.indA[,i])))
    }
  }
  # return(round(y[,dim(data.probA)[2]],6))
  return(y[,dim(data.probA)[2]])
}

#-----------------------------------------------------------------------------
# (DEPRECATED)
# return entire network matrix from indiv. covar (Var) + covariate itself as first column
# NetInd_k is a matrix (N x k) of network friend indicies; 
# When obs i doesn't have j friend, NetInd[i, j] = NA
#-----------------------------------------------------------------------------
.f.allCovars <- function(k, NetInd_k, Var, VarNm, misval = gvars$misXreplace) {
# .f.allCovars <- function(k, NetInd_k, Var, VarNm, misval = 0L) {
  assertthat::assert_that(is.matrix(NetInd_k))
  n <- length(Var)
  netVar_names <- NULL
  netVar_full <- NULL
  d <- matrix(0L, nrow = n, ncol = k + 1)
  d[, 1] <- Var
  d[, c(2 : (k+1))] <- apply(NetInd_k, 2, function(k_indx) {
                                          netVar <- Var[k_indx]  # netVar values for non-existing friends are automatically NA
                                          netVar[is.na(netVar)] <- misval
                                          return(netVar)
                                          })
  if (k>0) netVar_names <- netvar(varnm = VarNm, fidx = c(1:k))
  Var_names <- c(VarNm, netVar_names)
  colnames(d) <- Var_names
  return(d)
}

# (DEPRECATED)
#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
#---------------------------------------------------------------------------------
# METHOD 1 (empirical distribution of \bar{h}=\sum{h_i}) - NO LONGER USED
# SEE OLDER git repo for this implementation
  # For given data, estimate g[A|cA] and calculate est. of h_i(c) for given value of c and each i. 
  # * Draw B samples from emp distribution of cY, for each i;
  # * Assume the same dimensionality for cY acorss all i;
  # * Assume W's are iid, use g_N(A|W) and keep N_i constant;
  # * Drawing from distribution of g_N or g* to get A's;
  # * Calculate h_bar from the joint distribution of all c_i^Y over all i;
  
# METHOD 2 (fit logit to each P(A_i|A's,W's)) 
# Predict P(A_1=a1|W1,..,Wk)*P(A_2=a2|A_1,W1,..,Wk)*...*P(A_k=a_k|A_1,...,A_k,W1,...,Wk)
pred.hbars.old <- function(new_data=NULL, fit_h_reg_obj, NetInd_k) {
   pred_h_fcn <- function(h_logit_sep_k, m.gAi_vec) {
      .f_pred_h_k <- function(k, sel_k_indx, m.gAi_vec) {
        k_sel <- sel_k_indx
        A_nms_arr <- colnames(indA)[c(1:(k+1))]
        indA <- indA[,c(1:(k+1)), drop=FALSE]

        if (is.null(hform)) {
          W_nms_arr <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
          W_nms_arr <- W_nms_arr[!is.na(W_nms_arr)]         
        } else {
          W_nms_arr <- all.vars(as.formula(hform))[-1]
        }
        probAi_vec <- NULL
        for (k_i in (1:(k+1))) {
          Covars_nms <- c(A_nms_arr[-c(k_i:(k+1))], W_nms_arr)
          deterministic <- (determ_cols[,k_i] == 1L)
          predAs_idx <- (!deterministic & k_sel)
          if (!fit_fastmtx) {
            #************************************************                       
            probAi <- .f.gen.probA_N(data.frame(cY_mtx), !(predAs_idx), m.gAi_vec[[k_i]])
            #************************************************           
          }
          else {
            newdat_mtx <- cY_mtx[predAs_idx, Covars_nms, drop=FALSE]
            #************************************************           
            probAi <- rep_len(0,n)
            #************************************************
            if (sum(predAs_idx)>0) probAi[predAs_idx] <- .f_predict_fast(m.gAi_vec[[k_i]], newdat_mtx)
          }
          probAi_vec <- cbind(probAi_vec, probAi)
        }
        # return the likelihood for the entire vector of Ai's, based on cY_mtx (cY)
        return(list(h_vec = .f.cumprod.matrix(indA, probAi_vec)[k_sel]))
      }
      if (h_logit_sep_k) {
        #----------------------------------------------------------------------
        # VERSION 1: predict logit separately for each k in the network and combine results 
        h_vec <- rep_len(0,n)
        for (nFriends in (0:k)) {
          k_sel <- (new_data[,node_l$nFnode]==nFriends)
          # nothing to predict if no one observed with this size ntwk
          if (sum(k_sel)!=0) {
            h_vec[k_sel] <- .f_pred_h_k(k=nFriends, sel_k_indx=k_sel, m.gAi_vec[[nFriends+1]])$h_vec
          }
        }
        return(list(h_vec=h_vec))
      } else {
        #---------------------------------------------------------------------
        # VERSION 2: predict logistic to entire dataset with nFriends as an extra covariate
        return(.f_pred_h_k(k=k, sel_k_indx=rep_len(TRUE,n), m.gAi_vec))
      }
    }
    # ASSIGN VARIABLE NAMES BASED ON ARGUMENT fit_h_reg_obj
    k <- fit_h_reg_obj$k
    h_user <- fit_h_reg_obj$h_user
    h_user_fcn <- fit_h_reg_obj$h_user_fcn
    lbound <- fit_h_reg_obj$lbound
    max_npwt <- fit_h_reg_obj$max_npwt
    h_logit_sep_k <- fit_h_reg_obj$h_logit_sep_k
    node_l <- fit_h_reg_obj$node_l
    gform <- fit_h_reg_obj$gform
    hform <- fit_h_reg_obj$hform
    fit_fastmtx <- fit_h_reg_obj$fit_fastmtx
    netW_namesl <- fit_h_reg_obj$netW_namesl
    netA_names <-  fit_h_reg_obj$netA_names
    determ_cols_Friend <- fit_h_reg_obj$determ_cols_Friend
    # select only baseline covariates (W's) that are included in the original gform:

    if (!is.null(new_data)) {
      determ.g_user <- new_data$determ.g
      determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
      determ_cols <- (determ_cols_user | determ_cols_Friend)
    }
    if (is.null(new_data)) {    # use original fitted data for prediction
      new_data <- fit_h_reg_obj$cY_mtx_fitted      
      determ_cols <- fit_h_reg_obj$determ_cols_fitted
    }
    n <- nrow(new_data)
    indA <- as.matrix(new_data[, netA_names])
    if (is.null(hform)) {
      W_nms <- unlist(netW_namesl)
    } else {
      W_nms <- all.vars(as.formula(hform))[-1]
    }
    if (!(node_l$nFnode%in%W_nms)) { W_nms <- c(W_nms, node_l$nFnode) }
    cY_mtx <- cbind(determ_cols, as.matrix(new_data[, c(node_l$nFnode, W_nms, netA_names)]))
    #---------------------------------------------------------------------
    # MAIN BODY OF THE FUNCTION
    if (h_user==FALSE) {
      P.hbar.c <- pred_h_fcn(h_logit_sep_k, fit_h_reg_obj$m.gAi_vec_g)
      P.hbar.star.c <- pred_h_fcn(h_logit_sep_k, fit_h_reg_obj$m.gAi_vec_gstar)
    }
    else {
      h_bars_user <- h_user_fcn(k, data, node_l, NetInd_k)
      P.hbar.c <- h_bars_user$P.hbar.c
      P.hbar.star.c <- h_bars_user$P.hbar.star.c
    }

    h_tilde <- P.hbar.star.c$h_vec / P.hbar.c$h_vec
    h_tilde[is.nan(h_tilde)] <- 0     # 0/0 detection
    h_tilde <- bound(h_tilde, c(0,1/lbound))
    df_h_bar_vals <- data.frame(cY.ID = 0, 
                                h.star_c = P.hbar.star.c$h_vec, 
                                h_c = P.hbar.c$h_vec, 
                                h = h_tilde
                                )
    return(list(df_h_bar_vals=df_h_bar_vals))
}


# (DEPRECATED)
# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars.old <- function(data, h_fit_params) {
    #---------------------------------------------------------------------------------
    # 1) return observed network data cY_mtx if is.null(f.g_name)
    # 2) pass only one netW, which can be separate for g_star and g_0
    # SAMPLE A LARGE DATASET of cY's for given functions g* or g0, of size p*N for some p
    # Get rid of the loop, by assigning .f.allCovars ALL AT ONCE
    # NOTE (OS 06/02/15): The only reason we need to pass netW and netW_full is because g_0 and g^*
    # are assumed to be based on the same sW! This shouldn't be the case, need to allow sW_g ad sW_gstar to be different
    #---------------------------------------------------------------------------------
    get_hfit_data <- function(cY_mtx, k, Anode, NetInd_k, netW, f.g_name, f.g_args, p, misval = 0L)  {
      samp_g_data <- function(df_sel) {
        resamp_A <- f.gen.A.star(k, df_sel, f.g_name, f.g_args)
        resamp_netA <- .f.allCovars(k, NetInd_k, resamp_A, Anode, misval = misval)
        fit.g_data <- cbind(netW, resamp_netA)
        return(fit.g_data)
      }
      if (is.null(f.g_name)) { return(cY_mtx) }
      n <- nrow(netW)
      df_samp_g <- samp_g_data(netW_full)  # big resampled matrix of c's (repeated data p times)
      fit.g_data_large <- matrix(nrow=(n*p), ncol=ncol(df_samp_g))
      colnames(fit.g_data_large) <- colnames(df_samp_g)
      fit.g_data_large[c(1:n),] <- as.matrix(df_samp_g)
      for (i in (2:p)) {
        fit.g_data_large[c(((i - 1) * n + 1):(n * i)), ] <- as.matrix(samp_g_data(netW_full))
      }
      return(fit.g_data_large)
    }
    # defining the vector of c^A's that needs evaluation under h(c) 
    .f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))
    #---------------------------------------------------------------------
    # Estimate h with logistic loss fcn (regression-based)
    # Estimate P(A_1|W1,..,Wk)*P(A_2|A_1,W1,..,Wk)*..*P(A_k|A_1,...,A_k,W1,...,Wk)
    # fitting each component P(A_i) as a logistic regression, for g_0 and g*
    # Does not rely on g_N - fitted model for g_0!!!!
    # Saves the estimated models for \bar{h}
    #---------------------------------------------------------------------  
    fit_h_fcn <- function(fit.g_data, p, h_logit_sep_k, netW_namesl, indA) {
      #----------------------------------------------------------------------
      # VERSION 1: Fit logit separately for each k in the network and combine together 
      # call logistic loss est. of h for each nFriends=0,..,k separately
      # defined k and k_sel values
      # .... REMOVED ...
      #---------------------------------------------------------------------
      # VERSION 2 for logistic loss h:
      # add extra smoothing: fit logistic to entire dataset with nFriends as an extra covariate
      # save the estimated model for \bar{h}
      k_sel <- rep_len(TRUE, n)
      #---------------------------------------------------------------------
      A_nms_arr <- colnames(indA)[c(1:(k+1))]
      indA <- indA[,c(1:(k+1)), drop=FALSE]
      if (is.null(hform)) {
        stop("NOT IMPLEMENTED, SUPPLY hform")
        # W_nms_arr <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
        # W_nms_arr <- W_nms_arr[!is.na(W_nms_arr)]
      } else {
        W_nms_arr <- all.vars(as.formula(hform))[-1]
      }
      probAi_vec <- NULL
      m.gAi_vec <- NULL
      for (k_i in (1:(k+1))) {  # factorization by nFriends (k) value
        A_i_nm <- A_nms_arr[k_i] # A variable we are predicting
        Covars_nms <- c(A_nms_arr[-c(k_i:(k+1))], W_nms_arr) # dependent covars
        deterministic <- (determ_cols[,k_i] == 1L)   # fit only to non-deterministic trt nodes
        fitAs_idx <- (!deterministic & k_sel)
        # USE glm.fit FUNCTION FOR FASTER FITTING (using design matrix format)
        Y_vals <- fit.g_data[rep.int(fitAs_idx, p), A_i_nm]
        X_mat <- as.matrix(fit.g_data[rep.int(fitAs_idx, p), Covars_nms, drop=FALSE], byrow=FALSE)
        m.gAi <- .f.est_binom_fast(X_mat, Y_vals)
        newdat_mtx <- cY_mtx[fitAs_idx, Covars_nms, drop=FALSE]   # original data to predict A

        # ************************************************
        probAi <- rep_len(0,n)
        # ************************************************

        if (sum(fitAs_idx)>0) probAi[fitAs_idx] <- .f_predict_fast(m.gAi, newdat_mtx)
        m.gAi_vec <- c(m.gAi_vec, list(m.gAi))
        probAi_vec <- cbind(probAi_vec, probAi)
      }
      # likelihood for the entire vector of Ai's, based on cY_mtx (cY)
      h_vec <- .f.cumprod.matrix(indA, probAi_vec)[k_sel]
      return(list(h_vec=h_vec, m.gAi_vec=m.gAi_vec))
    }

    #---------------------------------------------------------------------------------
    # PARAMETERS FOR LOGISTIC ESTIMATION OF h
    #---------------------------------------------------------------------------------
    # replace with p adaptive to k: p <- 100*(2^k)
    n <- nrow(data)
    n_samp_g0gstar <- h_fit_params$n_samp_g0gstar 
    family <- h_fit_params$family
    k <- h_fit_params$Kmax
    node_l <- h_fit_params$node_l
    NetInd_k <- h_fit_params$NetInd_k
    lbound <- h_fit_params$lbound
    max_npwt <- h_fit_params$max_npwt

    h_logit_sep_k=h_fit_params$h_logit_sep_k
    h_user=h_fit_params$h_user; h_user_fcn=h_fit_params$h_user_fcn; hform=h_fit_params$hform
    # m.gN=h_fit_params$m.g0N # not used
    f.g.star=h_fit_params$f.g.star; f.g_args=h_fit_params$f.g_args
    f.g0=h_fit_params$f.g0; f.g0_args=h_fit_params$f.g0_args
    fit_fastmtx <- TRUE
    gform <- h_fit_params$gform
    #---------------------------------------------------------------------------------

    # select only baseline covariates (W and W_net) that are included in hform
    nFnode <- node_l$nFnode
    Anode <- node_l$Anode

    # W's aren't resampled (assumed independent but not iid) => get all network W's from the original data
    netW <- NULL
    netW_namesl <- list()
    hform_covars <- all.vars(as.formula(hform))[-1]
    for (Wnode in node_l$Wnodes) {
      netW_names <- NULL
      # New condition (02/16/14):
      # only add netW_i covars for bsl covariate W if W doesn't already have "net" in its name
      netWnm_true <- (length(agrep("netF", Wnode, max.distance=list(insertions=0, substitutions=0)))==1)
      if (k>0 & !netWnm_true) netW_names <- netvar(Wnode, c(1:k))
      allW_names <- c(Wnode, netW_names)
      # only include bsl covariates (W, netW) if model for g_N depends on any of them
      if (!netWnm_true) { # only get netW's for non-network covars:
        netW_i <- .f.allCovars(k, NetInd_k, data[, Wnode], Wnode)
        netW <- cbind(netW, netW_i)
        if (any(allW_names %in% hform_covars)) netW_namesl <- c(netW_namesl, list(colnames(netW_i)))
      } else {
        netW <- cbind(netW, as.matrix(subset(data, select=Wnode)))
        if (any(allW_names %in% hform_covars)) netW_namesl <- c(netW_namesl, list(Wnode))
      }
    }

    # OS 06/01/15: moved from fit_h_fcn():
    if (is.null(hform)) {
      W_nms <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
      W_nms <- W_nms[!is.na(W_nms)]
    } else {
      W_nms <- all.vars(as.formula(hform))[-1]
    }
    if (!(nFnode%in%W_nms)) { W_nms <- c(W_nms,nFnode) }

    # 02/10/14: was a bug: get all W's for generating g_star (otherwise can get an error)
    netW_full <- NULL
    for (Wnode in node_l$Wnodes) {
        netW_i <- .f.allCovars(k, NetInd_k, data[, Wnode], Wnode)
        netW_full <- cbind(netW_full, netW_i)
    }
    netW_full <- cbind(netW_full, as.matrix(subset(data, select = nFnode)))

    #-----------------------------------------------------------
    # ADDING DETERMINISTIC FLAG COLUMNS TO netW
    #-----------------------------------------------------------
    # get deterministic As for the entire network of each unit (set by user)
    determ.g_user <- data$determ.g
    # determ.gvals_user <- data[,Anode] # add the values for deterministic nodes as well
    determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
    determ_cols_Friend <- 1L - .f.allCovars(k, NetInd_k, rep_len(1L, n), "determ.g_Friend")
    determ_cols <- (determ_cols_user | determ_cols_Friend)
    print("head(determ_cols_Friend)"); print(head(determ_cols_Friend))
    print("head(determ_cols)"); print(head(determ_cols))

    #-----------------------------------------------------------
    # DEFINING LOGISTIC REGs AND THE SUBSETING EXPRESSIONS
    # (1 expression per 1 regression P(sA[j]|sA[j-1:0], sW))
    #-----------------------------------------------------------
    A_nms <- netvar(Anode, c(0:k))
    regs_idx <- seq_along(A_nms)-1
    indA <- .f.allCovars(k = k, NetInd_k = NetInd_k, Var = data[,Anode], VarNm = Anode, misval = 0L)
    netW <- cbind(determ_cols, netW, as.matrix(subset(data, select=nFnode)))
    cY_mtx <- cbind(netW, indA)
    #-----------------------------------------------------------
    # BODY OF MAIN FUNCTION:
    #-----------------------------------------------------------
    ##########################################
    message("fitting h under g_0...")
    ##########################################
    p_h0 <- ifelse(is.null(f.g0), 1, n_samp_g0gstar)
    fit.g0_dat <- get_hfit_data(cY_mtx = cY_mtx, k = k, Anode = Anode, NetInd_k = NetInd_k,
                                netW = netW, f.g_name = f.g0, f.g_args = f.g0_args, p = p_h0)
    fit.g0_dat <- data.frame(fit.g0_dat)
    P.hbar.c <- fit_h_fcn(fit.g_data = fit.g0_dat, p = p_h0, h_logit_sep_k = h_logit_sep_k, netW_namesl = netW_namesl, indA = indA)
    ##########################################
    message("fitting h under g_star...")
    ##########################################
    # produces a matrix (not df, so needs to be converted for faster glm.fit)
    t1 <- system.time(
      fit.gstar_dat <- get_hfit_data(k = k, Anode = Anode, NetInd_k = NetInd_k, 
                                      netW = netW, f.g_name = f.g.star, f.g_args = f.g_args,
                                      p = n_samp_g0gstar)
    )
    # determ_cols may change under gstar, e.g., when P_g^*(A=1|W)=1 for some W
    fit.gstar_dat <- data.frame(fit.gstar_dat)
    P.hbar.star.c <- fit_h_fcn(fit.g_data = fit.gstar_dat, p = n_samp_g0gstar, h_logit_sep_k = h_logit_sep_k, netW_namesl = netW_namesl, indA = indA)
    # ##########################################
    # 3) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
    ##########################################
    h_tilde <- P.hbar.star.c$h_vec / P.hbar.c$h_vec
    h_tilde[is.nan(h_tilde)] <- 0     # 0/0 detection
    h_tilde <- bound(h_tilde, c(0,1/lbound))
    df_h_bar_vals <- data.frame(cY.ID = .f.mkstrNet(cY_mtx),
                                h.star_c = P.hbar.star.c$h_vec,
                                h_c = P.hbar.c$h_vec,
                                h = h_tilde
                                )

    fit_h_reg_obj <- list(m.gAi_vec_g = P.hbar.c$m.gAi_vec,
                          m.gAi_vec_gstar = P.hbar.star.c$m.gAi_vec,
                          k = k,
                          lbound = lbound,
                          h_logit_sep_k = h_logit_sep_k,
                          h_user = h_user,
                          h_user_fcn = h_user_fcn,
                          fit_fastmtx = fit_fastmtx,
                          node_l = node_l,
                          netW_namesl=netW_namesl,
                          gform = gform,
                          hform = hform,
                          netA_names=colnames(indA),
                          determ_cols_Friend = determ_cols_Friend,
                          determ_cols_fitted = determ_cols,
                          cY_mtx_fitted = cY_mtx)
    return(list(df_h_bar_vals = df_h_bar_vals, fit_h_reg_obj = fit_h_reg_obj))
}


# (DEPRECATED)
#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
#---------------------------------------------------------------------------------
get_all_ests.old <- function(data, est_obj) {
  # n <- nrow(data)
  # node_l <- data$nodes
  node_l <- est_obj$node_l
  # nFnode <- node_l$nFnode
  # Anode <- node_l$Anode
	Ynode <- node_l$Ynode

	Y <- data[, Ynode]
	determ.Q <- data[, "determ.Q"]
  # determ.g <- data[, "determ.g"]
  QY.init <- data[, "QY.init"] # initial Q fit
  off <- qlogis(QY.init)  # offset

  # new way of getting params (when data is an DatNet.sWsA object)
  # node_l <- data$nodes

  #************************************************
  # ESTIMATORS
  #************************************************

  #************************************************
  # IPTW_h estimator (based on h^*/h_N clever covariate):
  #************************************************
  fit.hbars_t <- system.time(h_bars <- fit.hbars.old(data = data, h_fit_params = est_obj)) # fit the clever covariat
  # fit.hbars_t.new <- system.time(h_bars.new <- fit.hbars.new(data = data, h_fit_params = est_obj)) # fit the clever covariat
  df_h_bar_vals <- h_bars$df_h_bar_vals
  fit_h_reg_obj <- h_bars$fit_h_reg_obj
  h_iptw <- df_h_bar_vals$h
  Y_IPTW_h <- Y
  Y_IPTW_h[!determ.Q] <- Y[!determ.Q] * h_iptw[!determ.Q]
  # print("IPW Est (h)"); print(mean(Y_IPTW_h))

  print("time to fit h_bars"); print(fit.hbars_t)
  print("h est"); print(head(df_h_bar_vals))
  # print("time to fit h_bars new"); print(fit.hbars_t.new)
  # print("h est new"); print(head(h_bars.new$df_h_bar_vals))

  #************************************************
  # TMLE A: estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.star_reg_A <- glm(Y ~ -1 + h_iptw + offset(off), data = data.frame(Y = data[, Ynode], off = off, h_iptw = h_iptw),
                                						subset = !determ.Q, family = est_obj$family, control = ctrl), GetWarningsToSuppress(TRUE))
	QY.star <- Y
	if (!is.na(coef(m.Q.star_reg_A))) QY.star <- plogis(off + coef(m.Q.star_reg_A) * h_iptw)

  #************************************************
  # TMLE B: estimate the TMLE update via weighted univariate ML (espsilon is intercept)
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.star_reg_B <- glm(Y ~ offset(off), data = data.frame(Y = data[, Ynode], off = off), weights = h_iptw,
                                            subset = !determ.Q, family = est_obj$family, control = ctrl), GetWarningsToSuppress(TRUE))
  QY.star_B <- Y
  if (!is.na(coef(m.Q.star_reg_B))) QY.star_B <- plogis(off + coef(m.Q.star_reg_B))

  #************************************************
  # IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
	# 02/16/13: IPTW estimator (Y_i * prod_{j \in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
	g_iptw <- iptw_est(k = est_obj$Kmax, data = data, node_l = node_l, m.gN = est_obj$m.g0N,
                      f.g.star = est_obj$f.g.star, f.g_args = est_obj$f.g_args, family = est_obj$family,
                      NetInd_k = est_obj$NetInd_k, lbound = est_obj$lbound, max_npwt = est_obj$max_npwt,
                      f.g0 = est_obj$f.g0, f.g0_args = est_obj$f.g0_args)
  Y_IPTW_net <- Y
  Y_IPTW_net[!determ.Q] <- Y[!determ.Q] * g_iptw[!determ.Q]

  #************************************************
  # IPTW-based clever covariate TMLE (based on FULL likelihood factorization), covariate based fluctuation
  #************************************************
	SuppressGivenWarnings(m.Q.star_iptw <- glm(Y ~ -1 + g_iptw + offset(off),
                                						data = data.frame(Y = data[, Ynode], off = off, g_iptw = g_iptw),
                                						subset = !determ.Q, family = est_obj$family, control = ctrl),
                                						GetWarningsToSuppress(TRUE))

  parsubmodel_fits <- rbind(coef(m.Q.star_reg_A), coef(m.Q.star_reg_B), coef(m.Q.star_iptw))
  rownames(parsubmodel_fits) <- c("epsilon (covariate)", "alpha (intercept)", "iptw epsilon (covariate)")
  print("parsubmodel_fits"); print(parsubmodel_fits)

  #************************************************
  # Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
  # TO DO: create a MC specific object with defined structure: models on Q's, data info, etc...
  # TO DO: create a MC evaluation that can work for any g^*
	#************************************************
  MC_fit_params <- append(est_obj,
                      list(m.Q.star_h_A = m.Q.star_reg_A,
                          m.Q.star_h_B = m.Q.star_reg_B,
                          m.Q.star_iptw = m.Q.star_iptw,
                          hstar = df_h_bar_vals$h.star_c,
                          hgN = df_h_bar_vals$h_c,
                          h_tilde = df_h_bar_vals$h))

  # run M.C. evaluation estimating psi under g^*:
  syst1 <- system.time(MCS_res <- get.MCS_ests(data = data,  MC_fit_params = MC_fit_params, fit_h_reg_obj = fit_h_reg_obj))
	print("time to run MCS: "); print(syst1);

  #************************************************
  # TO DO: come up with a better way to handle various estimators below:
  psi_mle <- MCS_res[names(MCS_res) == "gcomp_mle"]
  psi_tmle_A <- MCS_res[names(MCS_res) == "tmle_A"]
  psi_tmle_B <- MCS_res[names(MCS_res) == "tmle_B"]
  psi_tmle_iptw <- MCS_res[names(MCS_res) == "tmle_iptw"]
  psi_iptw_h <- mean(Y_IPTW_h)  # IPTW estimator based on h - clever covariate
  psi_iptw <- mean(Y_IPTW_net)  # IPTW estimator based on full g factorization (prod(g))

  # COMPONENTS OF ASYMPTOTIC VARIANCE FOR TMLE_NET (E_g^*[\bar{Q}_0(A^s,W^s|W^s)]-\psi_0):
  # SUBSTRACT overall estimate of psi_0 from fW_i i-specific components
  fWi_init <- MCS_res[agrep("fWi_init_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_init_A <- fWi_init - psi_tmle_A
  fWi_init_B <- fWi_init - psi_tmle_B
  fWi_init_tmleiptw <- fWi_init - psi_tmle_iptw

  fWi_star_A <- MCS_res[agrep("fWi_star_A_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_star_B <- MCS_res[agrep("fWi_star_B_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_star_A <- fWi_star_A - psi_tmle_A
  fWi_star_B <- fWi_star_B - psi_tmle_B
  print("fWi_star_A and fWi_star_B"); print(c(fWi_star_A=mean(fWi_star_A), fWi_star_B=mean(fWi_star_B)));

  return(list( tmle_A = psi_tmle_A,
               tmle_B = psi_tmle_B,
               tmle_iptw = psi_tmle_iptw,
               iptw_h = psi_iptw_h,
               iptw = psi_iptw,
               mle = psi_mle,
               fWi_init_A = fWi_init_A, fWi_star_A = fWi_star_A,
               fWi_init_B = fWi_init_B, fWi_star_B = fWi_star_B,
               fWi_init_tmleiptw = fWi_init_tmleiptw,
               h_iptw = h_iptw, iptw_reg = g_iptw,
               QY.init = QY.init, QY.star = QY.star))
}
