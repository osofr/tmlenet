#' @title tmlenet-package
#' @docType package
#' @import R6

# data.table
#' @author Oleg Sofrygin, Mark J. van der Laan
#' @description ...TO BE COMPLETED...
#' @name tmlenet-package
NULL

#######################################################################
######### BETA VERSION - NOT FOR DISTRIBUTION #########################
#######################################################################

#------------------------------------
# v0.2 Adding R6 classes for h estimation and syntax for arbitrary summary measures
# 07/11/15
# NETWORK TMLE
# authors: Oleg Sofrygin <sofrygin@berkeley.edu> and Mark van der Laan <laan@berkeley.edu>
#------------------------------------
# Arguments:
  #   data - Dataset of data.frame class with named nodes, wide format
  #   Anode - Name of the treatment variable
  #   Wnodes - Name(s) of confounding covariates
  #   Ynode -  Name of the outcome variable
  #   nFnode - Variable for network size of each observation
  #   Kmax - Constant for maximum number of friends any observation can have (cannot be exceeded)
  #   IDnode - Optional subject identifier variable in data, if not supplied the network in NETIDnode is assumed to be indexed by row numbers
  #   NETIDnode - String variable in data identifying subject's network by IDs or row numbers (space separated)
  #   Qform -  Regression formula for outcome Y; when NULL (default) Y is regressed on all summary measures defined in sW and sA
  #   gform -  Regression formula for treatment mechanism, g
  #   hform - Regression formula for the clever covariate, P(A^s | W^s) / P(A^{*s} | W^s) used for estimating under g_0 and g^*.
  #   AnodeDET - Column name for indicators of deterministic values of A, coded as (TRUE/FALSE) or (1/0) (to be removed: TRUE sets A=0), (to add: the observations with AnodeDET=TRUE/1 are assumed to have constant value for their Anode)
  # (TO ADD)  AnodeDETfun - function that evaluates to TRUE for observations with deterministically assigned A (alternative to AnodeDET)
  #   YnodeDET - Column name for indicators of deterministic values of Y, coded as (TRUE/FALSE) or (1/0) (to be removed: TRUE sets Y=1), (to add: the observations with YnodeDET=TRUE/1 are assumed to have constant value for their Ynode)
  # (TO ADD)  YnodeDETfun - function that evaluates to TRUE for observations with deterministically assigned Y (alternative to YnodeDET)
  #   Q.SL.library - SuperLearner libraries for outcome, Q (NOT IMPLEMENTED)
  #   g.SL.library - SuperLearner libraries for treatment mechanism, g (NOT IMPLEMENTED)
  #   gbound - One value for symmetrical bounds on g(A|W), or a vector containing upper and lower bounds
  #   f_gstar1 - Function for dynamic treatment regimen, can take any variables in data as arguments
  #   args_f_g1star - Additional arguments that will be passed to f_gstar1
  #   f_gstar2 - Optional dynamic treatment regimen for contrasts
  #   args_f_g2star - Optional additional arguments to be passed to f_gstar2
  #   f_g0 - When known, a function for generating A under true treatment mechanism, g0. Used only during estimation of the clevel covariate, h, via logistic regression after sampling A from this function
  #   h_f.g0_args - Additional arguments to be passed to f_g0
  #   h_user - Should a user supplied function be used to calculate the clever covariate, h
  #   h_user_fcn - User supplied function to calculate the clever covariate, h
  #   h_logit_sep_k - Flag for fitting a separate logistic regression for each strata of nFnode, used during estimation of the clever covariate, h
  #   family - Family specification for regression models, defaults to binomial. CURRENTLY ONLY BINOMIAL FAMILY IS IMPLEMENTED.
  #   alpha - alpha-level for CI calculation
  #   n_MCsims = ceiling(sqrt(nrow(data))) - number of Monte-Carlo simulation to perform for evaluation of psi under gstar & h_gstar
  #   verbose - Flag for controlling printing of messages (NOT IMPLEMENTED)
#---------------------------------------------------------------------------------
# Network Specification:
  #   1.
  #   IDnode - optional subject identifier column in d
  #   NETIDnode - network string of ID's column in d (space separated)
  #   2.
  #   NETIDnodes - variables indicating network IDs, k columns in data with "NA" for none (NOT IMPLEMENTED)
#---------------------------------------------------------------------------------
# Value:
  # EY_gstar1 - population mean under g1.star
  # EY_gstar2 - population mean under g2.star (if g2.star is provided)
  # ATE - additive treatment effect for EY_gstar1-EY_g0 (default) or EY_gstar1-EY_gstar2 (if f_gstar2 is provided)
  # Each element in the estimates of these is itself a list containing
    # • psi - parameter estimate for network effect
    # • pvalue_tmle - two-sided p-value
    # • var_tmle_Q.init - Influence-curve based variance of TMLE estimate
    # • c_var_tmle_Q.star - Conservative variance of TMLE estimate
    # • cabs_var_tmle_Q.star - Even more conservative variance of TMLE estimate
    # • CI_tmle_Q.init - corresponding alpha-based confidence interval
    # • c_CI_tmle_Q.star - corresponding alpha-based confidence interval
    # • cabs_CI_tmle_Q.star - corresponding alpha-based confidence interval
    # • psi_tmle_iptw - parameter estimate under inefficient TMLE
    # • psi_iptw - parameter estimate under IPTW
    # • psi_mle - parameter estimate under MLE
#---------------------------------------------------------------------------------

#------------------------------------
# REVISIONS (06/15/2013)
#------------------------------------
  # 1) Changed code to accept any number of W's 
  # 2) Defined deterministic nodes for A and for Q
  # 3) Passing g*() (stochastic intervention) as an argument

#------------------------------------
# REVISIONS (08/05/2013)
#------------------------------------
  # 1) Logit h_bar calculation is now based on SAMPLED A from g* or g_N
  # 2) Logit h_bar W's are based on empirical (as observed, no sampling)
  # 3) Added new variable for # of times to sample from g* and g_N (h.iter.logis)
  # 4) Default h.iter.logis set to 100
  # 5) IC-BASED variance for tmle is calculated
  # 6) Logit h_bar is now based on (c) only, doesn't depend on g_N, 
      # the model is saved for the MC evalution part (much faster run time)

#------------------------------------
# REVISIONS (11/18/2013)
#------------------------------------
  # *) Fixed error with IPTW
  # *) Much faster asymptotic variance calculation
  # *) Calculate ATE and IC.VAR for ATE if passed g.2.star function in addition to g.1.star
  # *) Conservative estimate of the Var_tmle (_c) that doesn't depend on correct specification of Q.0
  # *) Ad-hoc conservative estimate of the Var_tmle (_cabs) that doesn't depend on the assumption of 
        # positive correlation between friends' outcomes
  # *) Can directly pass the function for g_0 (generates A's under g_0 if g_0 known)
  # *) Can directly pass the function for h/h* given c

#------------------------------------
# REVISIONS (11/24/2013)
#------------------------------------
  # *) Added matrix-based computation of logistic reg for h (much faster)
  # *) Added calculation of f_Wi based on QY.star (in addition to f_Wi based on QY.init)
  # *) Fixed error when fitting h with several covariates (W1, W2)
  # *) Fixed issues with dimensionality when predicting h, with X_mat only 1 observation
  # *) Removed calculation of D*_Wi since we don't use it for as.var anymore

#------------------------------------
# REVISIONS (01/02/2014)
#------------------------------------
  # *) Changed the fit formulas for h_bar to be the same as the formulas for g_N, i.e. E[A|W] is fit to the same
  #     regression formula for IPTW, TMLE_IPTW and TMLE_EFF.
  #     Before h was calculated with A regressing on ALL network baseline covariates

#------------------------------------
# REVISIONS (02/11/2014)
#------------------------------------
  # *) Fixed possible bug:
  #     When fitting h_0 the model was using ALL observataions to fit E[netA_1|A,W1,netW]
  #     This was leading to TMLE bias for some simulations when Q misspecified
  #     e.g., observations with k=0 have all netA_i=0 for all i; k=1 have netA_i=0 for i>1 and so on...
  #     those observations where still being used when fitting the model for netA_1|A,W, netA_2|netA_1,A,W and so on...
  #     This could lead to scewed weights, fixed to fit only on non-determ A's
  # *) Fixed a bug in TMLE-efficient which resulted from previous revision (01/02/2014)
  #     f.gstar() was no longer passed all of baseline covariates when estimating h
  #     As a result the sampled A under g.star were incorrect, leading to incorrect estimates of weights h

#------------------------------------
# REVISIONS (21/03/2014)
#------------------------------------
  # *) New argument iidW_flag, when set to FALSE no resampling of W in MC evaluation (estimates psi that is conditional on W)
  # *) Passing object of parameters to MC sim function
  # *) Added inference based on iid ICs for multiple point treatments
  #     ***IPTW and TMLE-IPTW have more terms in the cross-prod***
  #     ***MLE inference based on IC for TMLE-efficient***
  #     ***Made sure fWi component is estimated based on m.Q.star (as it includes psi_i compenent over i)
  # *) Simulations: saving CI width in addition to coverage 
  #     To be able to compare different ests of as var (for TMLE-eff)
  
  # *) Changed weight truncation to be set to max 5% of the total weights - STILL WORKING IT OUT
  # *) Simulations: run iid simulation with multiple point treatments
  #     To show this TMLE works on iid data and gives correct coverage (no cross products for var ests)
########################################################################################

#------------------------------------
# REVISIONS (11/05/2014)
#------------------------------------
  # *) New iptw estimator that is based on h.star/h_0 alone
  # *) New intercept=based TMLE
  # *) Added separate regression specification for \bar{h} and \bar{h}^*

#------------------------------------
# TO DO LIST:
#------------------------------------
  # *) (11/05/14) Make as. var estimates to be at least \sum{(D^*^2)(O_i)}, i.e. lower bounded by as. var in iid case,
  #     since there is no reason to expect to have more information in dep. data than we would in iid case?
  # *) Add another ad-hoc as var est where negative cros prod terms are set to 0: 
  #     - Less conservative compared to as var est with sum of abs cros prod terms
  # *) rewrite predict/fit functions for h_bar estimation to allow for SL fitting
  # *) add debug on/off switch variable for printing
  # *) Allow for data-adaptive weight truncation wrt minimization of MSE (taking max of 5% weight of total as truth)
  # *) Evaluate coverage for data-adaptive psi_n,0 (MC eval under g.star for each sample of W's)
  # *) Implement EE-based estimator (see paper)
  # *) Allow SL to fit Q0_N, g0_N and h0_bar (i.e P(A_j|A_{1},..,A_{j-1}, W_i\inF_j))

`%+%` <- function(a, b) paste0(a, b)

.onAttach <- function(...) {
  packageStartupMessage("tmlenet")
  packageStartupMessage("The tmlenet package is still in beta testing. Please do not distribute. Interpret results with caution.")
}

#-----------------------------------------------------------------------------
#General utilities
#-----------------------------------------------------------------------------

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred", "prediction from a rank-deficient fit may be misleading")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

CheckInputs <- function(data, inputparams) {    # Error checking for inputs
  # ....
}

# THIS NEEDS MORE EVALUATION, DOESN'T SEEM TO WORK AS INTENDED DURING MC EVALUTION
# (GET HUGE ABSOLUTE VALUE WEIGHTS, THAT HAVE tiny % CONTRIBUTION)
scale_bound <- function(weights, max_npwt, n) {
  # scale weights
  weight_prop <- (weights/sum(weights))   # % contribution to the total weight, scaled by n
  weight_prop_byn <- (weights/sum(weights))*n   # % contribution to the total weight, scaled by n
  # print("weight summary before trunc"); print(summary(weights))
  # print("weight_prop before trunc, %"); print(summary(weight_prop))
  # print("weight_prop before trunc, scaled by n"); print(summary(weight_prop_byn))

  while (any(weight_prop_byn > (max_npwt+5))) {
    weights[which(weight_prop_byn >= max_npwt)] <- (max_npwt / n) * sum(weights)
    weight_prop_byn <- (weights/sum(weights)) * n   # % contribution to the total weight, sclaed by n
    # print("weight summary after trunc"); print(summary(weights))
    # print("weight_prop after trunc, %"); print(summary(weights/sum(weights)))
    # print("weight_prop after trunc, scaled by n"); print(summary(weight_prop_byn))
  }
  return(weights)
}

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}

# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
# data is an DatNet.sWsA object
# #todo 22 (get_all_ests) +0: rename get_all_ests into fit.param.fluct
# #todo 23 (get_all_ests) +0: move MC estimation (and h estimation?) outside get_all_ests
#---------------------------------------------------------------------------------
get_all_ests <- function(datNetObs, est_params_list) {
  # datNetObs$det.Y             # TRUE/FALSE for deterministic Y's
  # datNetObs$noNA.Ynodevals    # actual observed Y's
  # m.Q.init$getoutvarnm        # reg outvar name (Ynode)
  # datNetObs$YnodeVals         # visible Y's with NA for det.Y
  # m.Q.init$getoutvarval       # Yvals used in prediction (with det.Y obs set to NA)
  # m.Q.init$getprobA1          # predictions (for non-DET Y)
  # m.Q.init$getsubset          # valid subset (!det.Y)
  # m.Q.init$reg                # regression class (Qreg)

  nodes <- datNetObs$nodes
  Y <- datNetObs$noNA.Ynodevals # actual observed Y's
  determ.Q <- datNetObs$det.Y
  m.Q.init <- est_params_list$m.Q.init

  QY.init <- datNetObs$noNA.Ynodevals # getting all node vals, inc. deterministic
  QY.init[!datNetObs$det.Y] <- m.Q.init$getprobA1[!datNetObs$det.Y] # getting predictions P(Y=1) for non-DET Y
  off <- qlogis(QY.init)  # offset

  #************************************************
  # ESTIMATORS:
  #************************************************

  #************************************************
  # h^*/h_N clever covariate:
  #************************************************
  fit.hbars_t <- system.time(fit.hbars.res <- fit.hbars(datNetObs = datNetObs, est_params_list = est_params_list)) # fit the clever covariat

  datNetObs <- fit.hbars.res$datNetObs
  DatNet.gstar <- fit.hbars.res$DatNet.gstar

  m.h.fit <- fit.hbars.res$m.h.fit

  dat_hest <- fit.hbars.res$dat_hest
  # h_gstar <- fit.hbars.res$dat_hest$h_gstar
  # h_gN <- fit.hbars.res$dat_hest$h_gN
  h_wts <- fit.hbars.res$dat_hest$h_gstar_gN

  print("time to fit new fit.hbars.res:"); print(fit.hbars_t)
  print("new h est:"); print(dim(dat_hest)); print(head(dat_hest))

  #************************************************
  # IPTW_h estimator:
  #************************************************  
  Y_IPTW_h <- Y
  Y_IPTW_h[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]

  #************************************************
  # TMLE A: estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET
  # #todo 19 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.starA
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.starA <- glm(Y ~ -1 + h_wts + offset(off), data = data.frame(Y = Y, off = off, h_wts = h_wts),
                                						subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))
	QY.star_A <- Y
	if (!is.na(coef(m.Q.starA))) QY.star_A <- plogis(off + coef(m.Q.starA) * h_wts)

  #************************************************
  # TMLE B: estimate the TMLE update via weighted univariate ML (espsilon is intercept)
  # #todo 20 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.starB
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.starB <- glm(Y ~ offset(off), data = data.frame(Y = Y, off = off), weights = h_wts,
                                            subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))
  QY.star_B <- Y
  if (!is.na(coef(m.Q.starB))) QY.star_B <- plogis(off + coef(m.Q.starB))

  #************************************************
  # IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
	# 02/16/13: IPTW estimator (Y_i * prod_{j \in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
	# g_wts <- iptw_est(k = est_params_list$Kmax, data = data, node_l = nodes, m.gN = est_params_list$m.g0N,
  #                      f.gstar = est_params_list$f.gstar, f.g_args = est_params_list$f.g_args, family = "binomial",
  #                      NetInd_k = est_params_list$NetInd_k, lbound = est_params_list$lbound, max_npwt = est_params_list$max_npwt,
  #                      f.g0 = est_params_list$f.g0, f.g0_args = est_params_list$f.g0_args)
  #  Y_IPTW_g <- Y
  #  Y_IPTW_g[!determ.Q] <- Y[!determ.Q] * g_wts[!determ.Q]
  Y_IPTW_g <- rep(0, length(Y))

  #************************************************
  # IPTW-based clever covariate TMLE (based on FULL likelihood factorization), covariate based fluctuation
  # #todo 21 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star_giptw
  #************************************************
	# SuppressGivenWarnings(m.Q.star_giptw <- glm(Y ~ -1 + g_wts + offset(off),
  #                                						data = data.frame(Y = Y, off = off, g_wts = g_wts),
  #                                						subset = !determ.Q, family = "quasibinomial", control = ctrl),
  #                                						GetWarningsToSuppress(TRUE))  
  parsubmodel_fits <- rbind(coef(m.Q.starA), coef(m.Q.starB)) # parsubmodel_fits <- rbind(coef(m.Q.starA), coef(m.Q.starB), coef(m.Q.star_giptw))
  rownames(parsubmodel_fits) <- c("epsilon (covariate)", "alpha (intercept)") # rownames(parsubmodel_fits) <- c("epsilon (covariate)", "alpha (intercept)", "iptw epsilon (covariate)")
  print("new parsubmodel_fits: "); print(parsubmodel_fits)

  #************************************************
  # Run Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
	#************************************************
  MC_fit_params <- append(est_params_list,
                      list(m.Q.starA = m.Q.starA,
                          m.Q.starB = m.Q.starB,
                          # m.Q.star_giptw = m.Q.star_giptw,
                          dat_hest = dat_hest
                          ))

  syst1 <- system.time(MCS_res <- get.MCS_ests(datNetObs = datNetObs, DatNet.gstar = DatNet.gstar, 
                                                MC_fit_params = MC_fit_params, m.h.fit = m.h.fit))
  print("time to run MCS: "); print(syst1);

  ests <- c(tmle_A = MCS_res[["tmle_A"]],
            tmle_B = MCS_res[["tmle_B"]],
            tmle_g_iptw = MCS_res[["tmle_g_iptw"]],
            h_iptw = mean(Y_IPTW_h),  # IPTW estimator based on h - clever covariate
            g_iptw = mean(Y_IPTW_g),  # IPTW estimator based on full g factorization (prod(g))
            mle = MCS_res[["mle"]])

  ests_mat <- matrix(0L, nrow = length(ests), ncol = 1)
  ests_mat[, 1] <- ests
  rownames(ests_mat) <- names(ests); colnames(ests_mat) <- "estimate"

  wts_mat <- matrix(0L, nrow = datNetObs$nobs, ncol = 2)
  colnames(wts_mat) <- c("h_wts", "g_wts")
  wts_mat[, "h_wts"] <- h_wts
  # wts_mat[, "g_wts"] <- g_wts

  # Components of asymptotic variance for tmle_net (E_g^*[\bar{Q}_0(A^s,W^s|W^s)]-\psi_0):
  # SUBSTRACT overall estimate of psi_0 from fW_i i-specific components
  fWi_mat <- matrix(0L, nrow = datNetObs$nobs, ncol = 3)
  colnames(fWi_mat) <- c("fWi_Qinit", "fWi_Qstar_A", "fWi_Qstar_B")
  fWi_mat[,"fWi_Qinit"] <- MCS_res[agrep("fWi_init_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_mat[,"fWi_Qstar_A"] <- MCS_res[agrep("fWi_star_A_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_mat[,"fWi_Qstar_B"] <- MCS_res[agrep("fWi_star_B_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]

  QY_mat <- matrix(0L, nrow = datNetObs$nobs, ncol = 3)
  colnames(QY_mat) <- c("QY.init", "QY.star_A", "QY.star_B")
  QY_mat[,] <- cbind(QY.init, QY.star_A, QY.star_B)

  print(c(fWi_init_A = mean(fWi_mat[,"fWi_Qinit"] - ests["tmle_A"]), 
        fWi_init_B = mean(fWi_mat[,"fWi_Qinit"] - ests["tmle_B"]),
        fWi_star_A = mean(fWi_mat[,"fWi_Qstar_A"] - ests["tmle_A"]),
        fWi_star_B = mean(fWi_mat[,"fWi_Qstar_B"] - ests["tmle_B"])
        ));

  print("new MC.ests vec: "); print(ests)
  print("new MC.ests mat: "); print(ests_mat)

  return(list( ests_mat = ests_mat,
               wts_mat = wts_mat,
               fWi_mat = fWi_mat,
               QY_mat = QY_mat
              ))
}

get_vars_fromlist <- function(varname, sVar.map) {
  if (varname %in% names(sVar.map)) {
    as.vector(sVar.map[[varname]])
  } else {
    varname
  }
}
# NEW INTERFACE FOR SPECIFYING hform, Qform, gform:
process_regform <- function(regform, sW.map = NULL, sA.map = NULL) {
  if (length(regform)==0L) {
    return(list(outvars =  as.vector(unlist(sA.map)), predvars = as.vector(unlist(sW.map))))
  } else {

    # Getting predictors (sW names):
    regformterms <- terms(regform)
    sW.names <- attributes(regformterms)$term.labels 
    sW.names.alt <- colnames(attributes(regformterms)$factors)
    assert_that(all(sW.names == sW.names.alt))

    # Getting outcomes (sA names):
    out.var <- rownames(attributes(regformterms)$factors)[1] # character string
    out.vars.form <- as.formula(". ~ " %+% out.var)
    out.vars.terms <- terms(out.vars.form)
    sA.names <- attributes(out.vars.terms)$term.labels

    outvars <- unlist(lapply(sA.names, get_vars_fromlist, sA.map))
    predvars <- unlist(lapply(sW.names, get_vars_fromlist, sW.map))
    return(list(outvars = outvars, predvars = predvars))
  }
}

# ------------------------------
# MAIN TMLE ESTIMATOR FUNCTION:
# ------------------------------
# ------------------------------
# 4 layers of model spec's:
# ------------------------------
# 1) spec only Anode & Wnodes => assumes sW = W[[0:Kmax]] & sA = A[[0:Kmax]]
# 2) spec sW & sA => assumes hform = "sA ~ sW", gform = "A ~ sW", Qform = "Y ~ sW + sA"
# 3) spec sW & sA and spec either of hform, gform & Qform => assumes hform.gstar = hform
# 4) 3 + spec separate hform.gstar (outcome have to be the same sA for both hform & hform.gstar)
# *) Note: sA & sW can be character vectors consisting of R expressions
tmlenet <- function(data, Kmax, Anode, AnodeDET = NULL, Wnodes, Ynode, YnodeDET = NULL,
                    nFnode = "nF", IDnode = NULL, NETIDnode,
                    f_gstar1, f_gstar2 = NULL,
                    sW = NULL, sA = NULL,

                    # Replacing Qform, gform & hform with sA,sW and hform.new, hform.gstar.new, Qform.new, gform.new
                    Qform = NULL, gform = NULL, hform = NULL,
                    Qform.new = NULL, hform.new = NULL, hform.gstar.new = NULL, gform.new = NULL,

                    verbose = FALSE, # NOT IMPLEMENTED

                    args_f_g1star = NULL, args_f_g2star = NULL, # REMOVE, NO NEED TO PASS ARGs to f_gstar1 or f_gstar2

                    opt.params = list( # "opt.params" arg list:
                      alpha = 0.05,
                      gbound = 0.005, 
                      family = "binomial", # NOT IMPLEMENTED
                      n_MCsims = ceiling(sqrt(nrow(data))),
                      onlyTMLE_B = TRUE,
                      f_g0 = NULL)
                    ) {

  gvars$verbose <- verbose

  #----------------------------------------------------------------------------------
  # ADDITIONAL ARGUMENTS (Removed from input args of tmlenet())
  #----------------------------------------------------------------------------------
  # onlyTMLE_B <- TRUE # if TRUE, only evalute the intercept TMLE (TMLE_B)
  iidW_flag <- FALSE
  Q.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction")
  g.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction")
  max_npwt <- 50
  h_logit_sep_k <- FALSE
  
  alpha <- ifelse(is.null(opt.params$alpha), 0.05, opt.params$alpha)
  gbound <- ifelse(is.null(opt.params$gbound), 0.005, opt.params$gbound)
  family <- ifelse(is.null(opt.params$family), "binomial", opt.params$family)
  n_MCsims <- ifelse(is.null(opt.params$n_MCsims), ceiling(sqrt(nrow(data))), opt.params$n_MCsims)
  onlyTMLE_B <- ifelse(is.null(opt.params$onlyTMLE_B), TRUE, opt.params$onlyTMLE_B)
  f_g0 <- NULL 
  if(!is.null(opt.params$f_g0)) f_g0 <- opt.params$f_g0

  #----------------------------------------------------------------------------------
  # MONTE-CARLO SIMULATION PARAMETERS
  #----------------------------------------------------------------------------------
  nQ.MCsims <- as.integer(n_MCsims)  # number of times to sample MC sim for Q (each of size n)
  ng.MCsims <- as.integer(n_MCsims)  # number of times to sample MC sim for h (each of size n)
  assert_that(is.count(nQ.MCsims))
  assert_that(is.count(ng.MCsims))

  max.err_est <- 0.1    # maximum percent error for MCS estimators

  #----------------------------------------------------------------------------------
  # PARAMETERS FOR ESTIMATING h under g0 & gstar
  #----------------------------------------------------------------------------------
  f.g0 <- f_g0
  f.g0_args <- NULL
  h_user_fcn <- NULL
  h_user <- !(is.null(h_user_fcn))
  #Q: 02/15/14: IS IT A GOOD IDEA TO HAVE THE SAME (OR ANY) UPPER BOUND ON g_star?
  if(length(gbound)==1) gbound <- c(gbound, 1 - gbound)

  #----------------------------------------------------------------------------------
  # Input checks:
  # #todo 8 (tmlenet) +0: check all names exist in data (Anode, Wnodes, Ynode, etc...)
  # #todo 11 (tmlenet) +0: Wnodes & Anode are no longer needed if provided sW, sA, give a warning that Wnodes,Anodes will be ignored
  #----------------------------------------------------------------------------------
  assert_that(is.data.frame(data))
  Kmax <- as.integer(Kmax)
  assert_that(is.count(Kmax))
  nobs <- nrow(data)
  node_l <- list(IDnode = IDnode, Anode = Anode, Wnodes = Wnodes, Ynode = Ynode, nFnode = nFnode, NETIDnode = NETIDnode)

  # Defining Deterministic Y and A node flags:
  if (is.null(AnodeDET)) {determ.g <- rep_len(FALSE, nobs)} else {determ.g <- (data[, AnodeDET] == 1)}
  if (is.null(YnodeDET)) {determ.Q <- rep_len(FALSE, nobs)} else {determ.Q <- (data[, YnodeDET] == 1)}

  # #todo 63 (tmlenet) +0: Replace all input data with this (to get rid of irrelevant columns in data)
  d_sel <- data.frame(subset(data, select = unlist(node_l)), determ.g = determ.g, determ.Q = determ.Q)

  #----------------------------------------------------------------------------------
  # Create an object with model estimates, data & network information that is passed on to estimation procedure
  #----------------------------------------------------------------------------------
  netind_cl <- NetIndClass$new(Odata = data, Kmax = Kmax, IDnode = IDnode, NETIDnode = NETIDnode, sep = ' ')
  NetInd_k  <- netind_cl$NetInd_k
  if (nFnode %in% names(data)) {
    message("performing nFnode consistency check: ")
    if (!all(as.integer(netind_cl$nF) == as.integer(data[,nFnode]))) {
      stop("column data[,nFnode] does not match automatically calculated values for nFnode.
            Either remove this column from data or re-calculate it correctly.")
    }
    message("...passed...")
  }

  #----------------------------------------------------------------------------------
  # Parse and evaluate the summary measures (in class Define_sVar):
  # #todo 17 (tmletnet) +5: pre-evaluate summary measures on small batch of data to get dims of sA & sW and to check for errors
    # Would need to sort out how to define NetInd_k for a subset of full data?
  # #todo 52 (tmlenet) +0: Accept sA & sW as character vectors / lists passed to tmlenet (in addition to current set-up)
    # When sW / sA are just lists of character vectors need to capture the calling env and call Define_sVar constructor:
      # user.env <- parent.frame()
      # user.env_l <- list(user.env = user.env)
      # sW <- do.call(Define_sVar$new, c(sW, list(type = "sW"), user.env_l))
      # sW.gstar <- do.call(Define_sVar$new, c(sW.gstar, list(type = "sW.gstar"), user.env_l))
      # sA <- do.call(Define_sVar$new, c(sA, list(type = "sA"), user.env_l))  
  # #todo 53 (tmlenet) +0: If no sVars were defined (default), use netW (Wnode[[0:Kmax]]) and netA for sVars (Anode[[0:Kmax]])
  # #todo 54 (tmlenet) +0: Check all sVar names are unique
  #----------------------------------------------------------------------------------

  # Testing the evaluation of summary measures:
  testm.sW <- sW$get.mat.sVar(data.df = data, netind_cl = netind_cl, addnFnode = node_l$nFnode)
  print("testm.sW"); print(head(testm.sW))
  print("testm.sW map"); print(sW$sVar.names.map)

  testm.sA <- sA$get.mat.sVar(data.df = data, netind_cl = netind_cl)
  print("testm.sA"); print(head(testm.sA))
  print("testm.sA map"); print(sA$sVar.names.map)

  #---------------------------------------------------------------------------------
  # BUILDING OBSERVED sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
  #---------------------------------------------------------------------------------
  datnetW <- DatNet$new(netind_cl = netind_cl, nodes = node_l, VarNodes = node_l$Wnodes, addnFnode = TRUE)
  datnetW$make.sVar(Odata = data, sVar.object = sW)
  print("head(obsdat.sW) as dat.sVar"); print(head(datnetW$dat.sVar))
  datnetW$fixmiss_sVar() # permanently replace NA values in sW with 0
  print("head(obsdat.sW) after fixmiss_sVar:"); print(head(datnetW$dat.sVar))
  datnetA <- DatNet$new(netind_cl = netind_cl, nodes = node_l, VarNodes = node_l$Anode)
  datnetA$make.sVar(Odata = data, sVar.object = sA)
  print("head(obsdat.sA) as dat.sVar"); print(head(datnetA$dat.sVar))
  # print("head(obsdat.sA)after fixmiss_sVar: "); print(head(datnetA$dat.sVar))
  # (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sA:
  message("cancelled adding DET nodes to sVar since all sVar are automatically get added to A ~ predictors + DETnodes...")
  # obsdat.sW <- O.datnetW$add_deterministic(Odata = data, userDETcol = "determ.g")$dat.sVar
  # print(head(obsdat.sW))
  obsYvals <- data[,node_l$Ynode]
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = obsYvals, det.Y = determ.Q)$make.dat.sWsA()
  print("head(datNetObs$dat.sWsA)"); print(head(datNetObs$dat.sWsA))

  # Testing NA for visible det.Y and true observed Y as protected:
  # browser()
  # stop()
  # determ.Q <- c(FALSE, FALSE, FALSE, rep(TRUE, length(determ.Q)-3))
  # length(determ.Q) == length(obsYvals)
  # datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = obsYvals, det.Y = determ.Q)$make.dat.sWsA()

  #----------------------------------------------------------------------------------
  # Optional regressions specs:
  # todo 58 (tmlenet, Q.sVars, g.sVars) +0: Check that outvars & predvars in Q.sVars & g.sVars actually exist in sW, sA
  #----------------------------------------------------------------------------------
  Q.sVars <- process_regform(as.formula(Qform.new), sW.map = c(sW$sVar.names.map, sA$sVar.names.map), sA.map = node_l$Ynode)
  print("Qform.new: " %+% Qform.new)
  print("Q.sVars"); print(str(Q.sVars))

  h.sVars <- process_regform(as.formula(hform.new), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
  print("hform.new: " %+% hform.new)
  print("h.sVars"); print(str(h.sVars))

  if (!is.null(hform.gstar.new)) {
    h.gstar.sVars <- process_regform(as.formula(hform.gstar.new), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
  } else {
    h.gstar.sVars <- h.sVars
  }
  print("hform.gstar.new: " %+% hform.gstar.new)
  print("h.gstar.sVars"); print(str(h.gstar.sVars))

  g.sVars <- process_regform(as.formula(gform.new), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
  print("gform.new: " %+% gform.new)
  print("g.sVars: "); print(str(g.sVars))
  #-----------------------------------------------------------
  # Defining and fitting regression for Y ~ sW + sA:
  #-----------------------------------------------------------
  # subset_expr <- try(parse(text="!misfun("%+%node_l$Ynode%+%")")[[1]])
  # #todo 45 (m.Q.init) +0: Allow the option of fitting a separate m.Q.init model for each nF value
  # print("m.Q.init predictors: "); print(Q.sVars$predvars)
  Qreg <- RegressionClass$new(outvar = node_l$Ynode,
                              predvars = Q.sVars$predvars,
                              subset = !determ.Q,
                              ReplMisVal0 = TRUE,
                              form = Qform)
  m.Q.init <- BinOutModel$new(glm = FALSE, reg = Qreg)$fit(data = datNetObs)$predict()

  # datNetObs$YnodeVals       # visible Y's with NA for det.Y
  # datNetObs$det.Y           # TRUE/FALSE for deterministic Y's
  # datNetObs$noNA.Ynodevals  # actual observed Y's
  # m.Q.init$getoutvarnm      # reg outvar name (Ynode)
  # m.Q.init$getoutvarval     # obsYvals after setting det.Y obs to NA
  # m.Q.init$getprobA1        # predictions (for non-DET Y)
  # m.Q.init$getsubset        # valid subset (!det.Y)
  # m.Q.init$reg              # regression class (Qreg)
  # print("summary(m.Q.init): "); print(summary(m.Q.init))

  # Not needed here, using only for comparison:
  # QY.init <- datNetObs$noNA.Ynodevals
  # QY.init[!datNetObs$det.Y] <- m.Q.init$getprobA1[!datNetObs$det.Y] # getting predictions P(Y=1) for non-DET Y

  #----------------------------------------------------------------------------------
  # When Qform is provided, use the formula based fit for Q.init instead of Qreg (sW+sA) fit
  #----------------------------------------------------------------------------------
  if (!is.null(Qform)) {
    net_d <- cbind(datNetObs$dat.sWsA, subset(data, select = node_l$Ynode))
    net_d[gvars$misfun(net_d)] <- gvars$misXreplace
    print("head(net_d)"); print(head(net_d, 5))

    m.Q.init.old <- f_est(net_d[!determ.Q,], Qform, family = family)
    QY.init.old <- data[, node_l$Ynode] # setting deterministic node values
    QY.init.old[!determ.Q] <- predict(m.Q.init.old, newdata = net_d[!determ.Q,], type = "response") # predict p(Y) for non determ nodes    
  
    print("new coef(m.Q.init): "); print(coef(m.Q.init))
    print("old coef(m.Q.init.old): "); print(coef(m.Q.init.old))
  }

  # dfcheck <- data.frame(QY.init = QY.init, QY.init.old = QY.init.old, diff = QY.init - QY.init.old)
  # head(dfcheck, 50)
  # browser()
  # stop()

  #----------------------------------------------------------------------------------
  # DEPRECATED... All regs are now defined via sW, sA summary measures...
  # Create net_d for fitting m.Q.init, m.g0N and m.h_g0, m.h_gstar
  # Check all RHS var names in Qform, gform, hform exist (NOT DONE)
  #----------------------------------------------------------------------------------
  # if (is.null(Qform)) {
  #   # default to main terms in sW and sA
  #   Qform <- node_l$Ynode %+% " ~ " %+%
  #             paste0(datNetObs$datnetW$names.sVar, collapse="+") %+%
  #             "+" %+%
  #             paste0(datNetObs$datnetA$names.sVar, collapse="+")
  # }
  # if (is.null(hform)) {
  #   hform <- "s.A ~ " %+% paste0(datNetObs$datnetW$names.sVar, collapse="+") # default to main terms in datNetObs$datnetW
  # }
  # TO BE REMOVED:
  if (is.null(gform)) {
    gform <- node_l$Anode %+% " ~ " %+% paste0(datNetObs$datnetW$names.sVar, collapse="+") # default to main terms in datNetObs$datnetW
  }

  message("Running tmlenet with... ");
  message("Qform: " %+% Qform)
  message("gform: " %+% gform)
  message("hform: " %+% hform)

  m.g0N.old <- f_est(net_d[!determ.g,], gform, family = family) # Set A = 0 when determ.g == 1
  print("coef(m.g0N.old)"); print(coef(m.g0N.old))
  d_sel <- cbind(d_sel, QY.init = QY.init.old) # (DEPRECATED, TO BE REMOVED)
  print("head(d_sel) old: "); print(head(d_sel))
	#----------------------------------------------------------------------------------
  # Defining and fitting regression for A ~ sW:
  # #todo 64 (tmlenet) +0: Fit the model for g_0(A,W) - move this to get_all_ests() to avoid confusion
	#----------------------------------------------------------------------------------
  greg <- RegressionClass$new(outvar = node_l$Anode,
                              predvars = g.sVars$predvars,
                              subset = !determ.g)
  m.g0N <- BinOutModel$new(glm = FALSE, reg = greg)$fit(data = datNetObs)
  # print("summary(m.g0N): "); print(summary(m.g0N))
  print("coef(m.g0N): "); print(coef(m.g0N))

  #----------------------------------------------------------------------------------
  # Create an object with model estimates, data & network information that is passed on to estimation procedure
  #----------------------------------------------------------------------------------
  # 1) define parameters for MC estimation of the substitution estimators
  # 2) define parameters for estimation of the efficient weights h(A^s|W^s)
  est_obj <- list(
                  onlyTMLE_B = onlyTMLE_B,
                  lbound = gbound[1],
                  max.err_eps = max.err_est,  # error tolerance for the mean/var M.C. estimate

                  m.g0N = m.g0N, m.Q.init = m.Q.init, # now BinOutModel objects

                  f.g0 = f.g0, 
                  sW = sW, sA = sA,

                  Q.sVars = Q.sVars, 
                  h.sVars = h.sVars, 
                  h.gstar.sVars = h.gstar.sVars,
                  g.sVars = g.sVars,

                  nQ.MCsims = nQ.MCsims, 
                  ng.MCsims = ng.MCsims,

                  h_logit_sep_k = h_logit_sep_k, # NOT IMPLEMENTED
                  h_user = h_user, # NOT IMPLEMENTED
                  h_user_fcn = h_user_fcn, # NOT IMPLEMENTED

                  max_npwt = max_npwt, # NOT IMPLEMENTED  # cap the prop weights scaled at max_npwt (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000)
                  Kmax = Kmax, # REMOVE, already saved in DatNet
                  node_l = node_l, #REMOVE, already saved in DatNet
                  NetInd_k = NetInd_k, #REMOVE, already saved in DatNet
                  family = family, #REMOVE, already saved in DatNet regs
                  f.g0_args = f.g0_args, #REMOVE, no longer used
                  gform = gform, #REMOVE, no longer used
                  hform = hform, #REMOVE, no longer used
                  iidW_flag = iidW_flag, #REMOVE (NOT USED)
                  n_MCsims = nQ.MCsims,         # REMOVE, keeping for now for compatibility
                  n_samp_g0gstar = ng.MCsims  # REMOVE, keeping for now for compatibility
                  )

  est_obj_g1 <- append(est_obj,
                      list(
                        f.gstar = f_gstar1,
                        f.g_args = args_f_g1star #REMOVE, no longer used
                        )
                      )

  if (!is.null(f_gstar2)) {
    est_obj_g2 <- append(est_obj,
                      list(
                        f.gstar = f_gstar2,
                        f.g_args = args_f_g2star #REMOVE, no longer used
                        )
                      )
  }

  #----------------------------------------------------------------------------------
	# DEPRECATED: Run TMLE univariate fluctuations for each g.star and/or ATE:
	#----------------------------------------------------------------------------------
  est_obj_g1$m.g0N <- m.g0N.old
  est_obj_g1$m.Q.init <- m.Q.init.old
	tmle_g1_out.old <- get_all_ests.old(data = d_sel, est_obj = est_obj_g1)
  
  tmle_g2_out.old <- NULL
  if (!is.null(f_gstar2)) {
    est_obj_g2$m.g0N <- m.g0N.old
    est_obj_g2$m.Q.init <- m.Q.init.old
    tmle_g2_out.old <- get_all_ests.old(data = d_sel, est_obj = est_obj_g2)
  }

  #----------------------------------------------------------------------------------
  # Run MC evaluation for substitution TMLE ests
  #----------------------------------------------------------------------------------
  est_obj_g1$m.g0N <- m.g0N
  est_obj_g1$m.Q.init <- m.Q.init
  tmle_g1_out <- get_all_ests(datNetObs = datNetObs, est_params_list = est_obj_g1)

  tmle_g2_out <- NULL
  if (!is.null(f_gstar2)) {
    est_obj_g2$m.g0N <- m.g0N
    est_obj_g2$m.Q.init <- m.Q.init
    tmle_g2_out <- get_all_ests(datNetObs = datNetObs, est_params_list = est_obj_g2)
  }

  #----------------------------------------------------------------------------------
  # Create output list (estimates, as. variances, CIs)
  #----------------------------------------------------------------------------------
  EY_gstar1 <- make_EYg_obj(alpha = alpha, onlyTMLE_B = onlyTMLE_B, datNetObs = datNetObs, tmle_g_out = tmle_g1_out)
  # EY_gstar1 <- make_EYg_obj(alpha = alpha, onlyTMLE_B = onlyTMLE_B, datNetObs = datNetObs, tmle_g_out = tmle_g1_out.old)

  EY_gstar2 <- NULL
  ATE <- NULL	
  if (!is.null(f_gstar2)) {
    EY_gstar2 <- make_EYg_obj(alpha = alpha, onlyTMLE_B = onlyTMLE_B, datNetObs = datNetObs, tmle_g_out=tmle_g2_out)
    ATE <- make_EYg_obj(alpha = alpha, onlyTMLE_B = onlyTMLE_B, datNetObs = datNetObs, tmle_g_out = tmle_g1_out, tmle_g2_out = tmle_g2_out)
	}

	tmlenet.res <- list(EY_gstar1 = EY_gstar1, EY_gstar2 = EY_gstar2, ATE = ATE)
	class(tmlenet.res) <- c(class(tmlenet.res), "tmlenet")

  print("Estimates w/ MC eval:"); print(tmlenet.res$EY_gstar1$estimates)
	return(tmlenet.res)
}