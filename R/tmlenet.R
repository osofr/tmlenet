#' @title tmlenet-package
#' @docType package
#' @import R6

# data.table
#' @author Oleg Sofrygin, Mark J. van der Laan
#' @description ...TO COMPLETE...
#' @name tmlenet-package
NULL

#######################################################################
######### BETA VERSION - NOT FOR DISTRIBUTION #########################
#######################################################################

#------------------------------------
# v0.2 Adding R6 classes for h estimation and syntax for arbitrary summary measures
# 05/18/14 
# NETWORK TMLE
# authors: Oleg Sofrygin <sofrygin@berkeley.edu> and Mark van der Laan <laan@berkeley.edu>
#------------------------------------
# Arguments:
  #   data - Dataset of data.frame class with named nodes, wide format
  #   Anode - Name of the treatment variable
  #   Wnodes - Name(s) of confounding covariates
  #   Ynode -  Name of the outcome variable
  #   nFnode - Variable for network size of each observation
  #   k_max - Constant for maximum number of friends any observation can have (cannot be exceeded)
  #   IDnode - Optional subject identifier variable in data, if not supplied the network in NETIDnode is assumed to be indexed by row numbers
  #   NETIDnode - String variable in data identifying subject's network by IDs or row numbers (space separated)
  #   Qform -  Regression formula for outcome, Q 
  #   gform -  Regression formula for treatment mechanism, g
  #   h_form - Regression formula for the clever covariate, P(A^s | W^s) / P(A^{*s} | W^s) used for estimating under g_0 and g^*.
  #   AnodeDET - Column name for indicators of deterministic values of A, coded as (TRUE/FALSE) or (1/0) (to be removed: TRUE sets A=0), (to add: the observations with AnodeDET=TRUE/1 are assumed to have constant value for their Anode)
  # (TO ADD)  AnodeDETfun - function that evaluates to TRUE for observations with deterministically assigned A (alternative to AnodeDET)
  #   YnodeDET - Column name for indicators of deterministic values of Y, coded as (TRUE/FALSE) or (1/0) (to be removed: TRUE sets Y=1), (to add: the observations with YnodeDET=TRUE/1 are assumed to have constant value for their Ynode)
  # (TO ADD)  YnodeDETfun - function that evaluates to TRUE for observations with deterministically assigned Y (alternative to YnodeDET)
  #   Q.SL.library - SuperLearner libraries for outcome, Q (NOT IMPLEMENTED)
  #   g.SL.library - SuperLearner libraries for treatment mechanism, g (NOT IMPLEMENTED)
  #   gbound - One value for symmetrical bounds on g(A|W), or a vector containing upper and lower bounds
  #   f.g1.star - Function for dynamic treatment regimen, can take any variables in data as arguments
  #   f.g1_args - Additional arguments to be passed to f.g1.star
  #   f.g2.star - Optional dynamic treatment regimen for contrasts
  #   f.g2_args - Optional additional arguments to be passed to f.g2.star
  #   h_f.g0 - When known, a function for generating A under true treatment mechanism, g0. Used only during estimation of the clevel covariate, h, via logistic regression after sampling A from this function
  #   h_f.g0_args - Additional arguments to be passed to h_f.g0
  #   h_user - Should a user supplied function be used to calculate the clever covariate, h
  #   h_user_fcn - User supplied function to calculate the clever covariate, h
  #   h_logit_sep_k - Flag for fitting a separate logistic regression for each strata of nFnode, used during estimation of the clever covariate, h
  #   W_indep - Flag for independence of confounding covariates, when set to FALSE new TMLE is estimated, conditional on all W's (NOT IMPLEMENTED)
  #   family - Family specification for regression models, defaults to binomial. CURRENTLY ONLY BINOMIAL FAMILY IS IMPLEMENTED.
  #   alpha - alpha-level for CI calculation
  #   verbose - Flag for controlling printing of messages (NOT IMPLEMENTED)
  #   n_MCsims=ceiling(sqrt(nrow(data))) - number of Monte-Carlo simulation to perform for evaluation of psi under g^*
  #   n_samp_g0gstar=20 - number of Monte-Carlo simulations to perform for evaluation of \h^* under g^*
#---------------------------------------------------------------------------------
# Network Specification:
  #   1.
  #   IDnode - optional subject identifier column in d
  #   NETIDnode - network string of ID's column in d (space separated)
  #   2.
  #   NETIDnodes - variables indicating network IDs, k columns in data with "NA" for none (NOT IMPLEMENTED)
#---------------------------------------------------------------------------------
# Value:
  # EY_g1.star - population mean under g1.star
  # EY_g2.star - population mean under g2.star (if g2.star is provided)
  # ATE - additive treatment effect for EY_g1.star-EY_g0 (default) or EY_g1.star-EY_g2.star (if f.g2.star is provided)
  # Each element in the estimates of these is itself a list containing
    # • psi_tmle - parameter estimate for network TMLE
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
#	........
# Installing:
# install.packages("/Users/monikaizano/Dropbox/Network_TMLE/TMLENET_package/Test_Package/tmlenet_0.1.tar.gz", repos = NULL, type="source")

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
  # *) rewrite predict/fit functions for h_bar estimation to accept any regression form and allow for SL fitting
  # *) add debug on/off switch variable for printing

  # *) Allow for data-adaptive weight truncation wrt minimization of MSE (taking max of 5% weight of total as truth)
  # *) Implement new TMLE with:
  #     - A^s and W^s (arbitrary summary functions of A's and W's in friends network)
  #     - Common g_i's accross all i's: g(A_i|W^s)
  #     - h_0=P(A^s|W^s), over all i's (fit as iid over pooled observations A^s,W^s)
  # *) Evaluate coverage for data-adaptive psi_n,0 (MC eval under g.star for each sample of W's)
  # *) Implement EE-based estimator (see paper)
  # *) Add influence curve-based variance calculation for EE-based estimator
  # *) Allow SL to fit Q0_N, g0_N and h0_bar (i.e P(A_j|A_{1},..,A_{j-1}, W_i\inF_j))
  # *) Explore different types of network formation mechanisms (Simulations)
  #     - Currently only looking at a completely randomly formed network, such design is UNREALISTIC

  # *) PROBLEM (09/29/14) with MC evaluation of psi under g_star:
    # *) when W's are resampled (iidW=TRUE) deterministic Y's can change (as a function of W)
    # *) possible solutions: a) only resample W for non-deterministic Y's
    # *) possible solutions: b) **** have a function that re-assigns deterministic Y's based on current vector of W's ****

`%+%` <- function(a, b) paste0(a, b)
# jpaste <- function(...) paste(..., sep="")

.onAttach <- function(...) {
  packageStartupMessage("tmlenet")
  packageStartupMessage("The tmlenet package is still in beta testing. Please do not distribute. Interpret results with caution.")
}

# Get summary measures for one or two tmlenet objects 
#(SEs, p-values, CIs)
# If two objects, include effect measures (additive effect, relative risk, odds ratio)
# summary.tmlenet <- function(object, control.object=NULL, 
						# estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  # #object is treatment, control.object is control
  # if (! is.null(control.object) && class(control.object) != "tmlenet") 
        # stop("the control.object argument to summary.tmlenet must be of class tmlenet")
  # if (! estimator %in% c("tmle", "iptw", "gcomp", "naive")) stop("estimator should be one of: tmle, iptw, gcomp, naive")
  # treatment.summary <- NULL
  # if (! is.null(control.object)) {
    # control.summary <- GetSummary(control.object$estimates[estimator], control.object$IC[[estimator]], loggedIC=FALSE)
    # effect.measures <- GetEffectMeasures(est0=control.object$estimates[estimator], 
                                          # IC0=control.object$IC[[estimator]], 
                                          # est1=object$estimates[estimator], 
                                          # IC1=object$IC[[estimator]])
    # effect.measures.summary <- lapply(effect.measures, function (x) GetSummary(x$est, x$IC, x$loggedIC))
  # } else {
    # control.summary <- effect.measures.summary <- NULL
  # }
  # ans <- list(treatment=treatment.summary, control=control.summary, 
              # effect.measures=effect.measures.summary, treatment.call=object$call, 
              # control.call=control.object$call, estimator=estimator)
  # class(ans) <- "summary.tmlenet"
  # return(ans)
# }


# Get summary measures for tmlenet parameters (standard errors, p-values, confidence intervals)
summary.tmlenet <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) 
{
  if (! estimator %in% c("tmle", "iptw", "gcomp")) stop("estimator should be one of: tmle, iptw, gcomp")  
  n <- nrow(IC)
  v <- apply(IC, 2, var)
  std.dev <- sqrt(v/n)
  pval <- 2 * pnorm(-abs(estimate / std.dev))
  CI <- GetCI(estimate, std.dev)
  cmat <- cbind(estimate, std.dev, CI, pval)
  dimnames(cmat) <- list(names(estimate), c("Estimate", "Std. Error", "CI 2.5%", "CI 97.5%", "p-value"))
  ans <- list(cmat=cmat, estimator=estimator)
  class(ans) <- "summary.tmlenet"
  return(ans)
}


# Print method for summary.tmlenet
print.summary.tmlenet <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  if (is.null(x$control)) {
    PrintCall(x$treatment.call)
    PrintSummary(x$treatment)
  } else {
    cat("\nControl ")
    PrintCall(x$control.call)
    cat("Treatment ")
    PrintCall(x$treatment.call)
    cat("Treatment Estimate:\n")
    PrintSummary(x$treatment)
    cat("\nControl Estimate:\n")
    PrintSummary(x$control)
    cat("\nAdditive Effect:\n")
    PrintSummary(x$effect.measure$ATE)
  }
  invisible(x)
}

# Print method for tmlenet
print.tmlenet <- function(x, ...) {
  PrintCall(x$call)
  cat("TMLE Estimate: ", x$estimates["tmle"], "\n")
  invisible(x)
}

# Print a call
PrintCall <- function(cl) {
  cat("Call:  ", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

# Print estimate, standard error, p-value, confidence interval
PrintSummary <- function(x) {
  cat("   Parameter Estimate: ", signif(x$estimate, 5))
  cat("\n    Estimated Std Err: ", signif(x$std.dev^2, 5))
  cat("\n              p-value: ", ifelse(x$pvalue <= 2*10^-16, "<2e-16",signif(x$pvalue, 5)))
  cat("\n    95% Conf Interval:",paste("(", signif(x$CI[1], 5), ", ", signif(x$CI[2], 5), ")", sep=""),"\n")
  invisible(x)
}

# Calculate estimate, SE, p-value, confidence interval
GetSummary <- function(estimate, IC, loggedIC) {
  if (is.null(IC)) {
    std.dev <- NA
  } else {
    n <- length(IC)
    std.dev <- sqrt(var(IC) / n)
  }
  if (loggedIC) {
    pvalue <- 2 * pnorm(-abs(log(estimate) / std.dev))
    CI <- exp(GetCI(log(estimate), std.dev))
  } else {
    pvalue <- 2 * pnorm(-abs(estimate / std.dev))
    CI <- GetCI(estimate, std.dev)
  }
  
  return(list(estimate=estimate, std.dev=std.dev, pvalue=pvalue, CI=CI))
}

# Calculate 95% confidence interval
GetCI <- function(estimate, std.dev) {
  x <- qnorm(0.975) * std.dev
  CI <- cbind("2.5%"=estimate - x, "97.5%"=estimate + x)
  return(CI)
}

# Calculate Average Treatment Effect
GetEffectMeasures <- function(est0, IC0, est1, IC1) {  
  names(est0) <- names(est1) <- NULL
  ATE <- est1 - est0
  if (is.null(IC0)) {
    ATE.IC <- NULL
  } else {
    ATE.IC <- -IC0 + IC1   
  }
  return(list(ATE=list(est=ATE, IC=ATE.IC, loggedIC=FALSE)))
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

}
# Run GLM or SuperLearner (in future vs) for fitting initial Q and g
# old call: (data, form, family)
Estimate <- function(form, d, subs, family, newdata, SL.library) {
  f <- as.formula(form)  
  if (any(is.na(d[subs, LhsVars(f)]))) stop("NA in Estimate")
  # stats <- CalcRegressionStats(d, f, subs)
  if (is.null(SL.library) || length(RhsVars(f)) == 0) { #in a formula like "Y ~ 1", call glm
    #estimate using GLM
    if (sum(subs) > 1) {
      SuppressGivenWarnings({
        m <- glm(f, data=d, subset=subs, family=family, control=glm.control(trace=FALSE, maxit=1000))
        predicted.values <- predict(m, newdata=newdata, type="response")
      }, GetWarningsToSuppress())
    } else {
      #glm breaks when sum(subs) == 1
      predicted.values <- rep(d[subs, LhsVars(f)], nrow(newdata))
    }
  } else {
    #estimate using SuperLearner
    if (family == "quasibinomial") family <- "binomial"
    #remove aliased columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
    rhs <- setdiff(RhsVars(f), rownames(alias(f, data=d[subs,])$Complete))  
    #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    new.subs <- apply(newdata[, rhs, drop=FALSE], 1, function (x) !any(is.na(x)))  
    
    m <- SuperLearner(Y=d[subs, LhsVars(f)], X=d[subs, rhs, drop=FALSE], SL.library=SL.library, 
                      verbose=FALSE, family=family, newX=newdata[new.subs, rhs, drop=FALSE])
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- m$SL.predict
  }
  return(list(values=predicted.values, model=m))
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

# NO LONGER USED IN NEW h PREDICTION METHOD
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

# NO LONGER USED IN NEW h PREDICTION METHOD
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

#---------------------------------------------------------------------------------
# G-Comp & TMLEs: Use Monte-Carlo to estimate psi under stochastic g^* 
#---------------------------------------------------------------------------------
  # For given data, take Q[Y|cY]=m.Q.init and calcualte est. of psi under g*=f.g.star using Monte-Carlo integration:
  # * W_i can be iid or not (in latter case W are not resampled);
  # * Draw from the distributions of W and g*(A|W), keeping N_i constant, recalculate cY and cA each time;
  # * Recalculate Y^c under g^*
  # * Repeat nrep times and average.
#---------------------------------------------------------------------------------
get.MCS_ests <- function(data, data_net, MC_fit_params, fit_h_reg_obj) {
  family <- MC_fit_params$family

  alphaTMLE_B_only <- MC_fit_params$alphaTMLE_B_only
  n <- nrow(data)
  k <- MC_fit_params$k  
  n_MCsims <- MC_fit_params$n_MCsims
  lbound <- MC_fit_params$lbound
  max_npwt <- MC_fit_params$max_npwt  
  max.err_eps <- MC_fit_params$max.err_eps

  node_l <- MC_fit_params$node_l
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
  Ynode <- node_l$Ynode
  iidW_flag <- MC_fit_params$iidW_flag
  NetInd_k <- MC_fit_params$NetInd_k

  m.Q.init <- MC_fit_params$m.Q.init
  m.Q.star_h_A <- MC_fit_params$m.Q.star_h_A
  m.Q.star_h_B <- MC_fit_params$m.Q.star_h_B
  m.Q.star_iptw <- MC_fit_params$m.Q.star_iptw
  m.g0N  <- MC_fit_params$m.g0N
  f.g.star <- MC_fit_params$f.g.star; f.g_args <- MC_fit_params$f.g_args
  f.g0 <- MC_fit_params$f.g0; f.g0_args <- MC_fit_params$f.g0_args
  # hstar <- MC_fit_params$hstar; hgN <- MC_fit_params$hgN; h_tilde <- MC_fit_params$h_tilde # not used

  # IID TMLE UPDATE
  m.iidQ.star_reg_B <- MC_fit_params$m.iidQ.star_reg_B
  m.iidQ.init <- MC_fit_params$m.iidQ.init

  # 02/07/12 Eliminated the need to make two passes through the data for monte-carlo # Creates one matrix of all Q(Y|W_i)
  .f_ests_all <- function(NetInd_k, emp_netW) {
  # .f_ests_all <- function(NetInd_k, NetInd_k_D_Wi, NetVec_D_fullsamp) {
    .f.gen.reps <- function(nrep, NetInd_k, emp_netW) 	{
    # .f.gen.reps <- function(nrep, NetInd_k, emp_netW, i_vec, fixW, NetReInd_k, nsubsamp, nF_vec, j_i_in_Fj, j_fullsamp)   {  
      # Generate full sample   c=(A,W) under g* (and iid or not for Ws)
       .f.gen.sample <- function(NetInd_k, n) { 
        if (iidW_flag) {  # Version 1: Resample W's with replacement
          resamp_idx <- sample(c(1:n), n, replace=TRUE)
          netW <- NULL # get all network W's from the original data under resampled IDs
          for (Wnode in node_l$Wnodes) {
            netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[resamp_idx,Wnode], Wnode)))
          }
          full_dfW <- netW
        } else {  # Version 2: No resampling of W's (W's are indep. but not iid, using NPMLE that puts mass 1 on all obs i=1,...,N)
          resamp_idx <- c(1:n)
          full_dfW <- emp_netW
        }
        nFriend <- subset(data, select=nFnode) # nFriends (network) is never resampled        
        # print("full_dfW"); print(head(full_dfW,10))
  			resamp_A <- f.gen.A.star(k, data.frame(full_dfW, subset(data, select=nFnode)), f.g.star, f.g_args)
  			full_dfA <- .f.allCovars(k, NetInd_k, resamp_A, Anode)
  			determ_df <-  data.frame(determ.g=data$determ.g[resamp_idx], determ.Q=data$determ.Q[resamp_idx]) # get deterministic nodes, also resampled, since fcn of W's
        Y_resamp <- subset(data, select=Ynode) # use observed Y's - INCORRECT, BASED ON NEW W (iidW=TRUE), deterministic Y's might change values
        
        resamp_data <- data.frame(full_dfW, full_dfA, Y_resamp, nFriend, determ_df)
        # print("resamp_data"); print(head(resamp_data))
  			return(resamp_data)
  		}

      # MLE - Predict E[Y|g_star] (QY.init) for each i, based on the initial model for E[Y|C^Y] (m.Q.init)
      .f_pred_barQ.init <- function(m.Q.init, samp_data) {
        #deterministic nodes for Q
        determ.Q <- samp_data$determ.Q
        # MLE Subs Estimator (est probY based on model for Q_Y)
        # predict only for those not infected at bsl, W2==0
        #************************************************
        # QY <- rep_len(1, nrow(samp_data))
        QY <- predict(m.Q.init, newdata=samp_data, type="response")
        QY[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect when W's are resampled
        #************************************************
        # QY[!determ.Q] <- predict(m.Q.init, newdata=samp_data[!determ.Q,], type="response")
        return(QY)
      }
      # TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the coefficient epsilon update model for E[Y|C^Y] (m.Q.star_h_A)      
      .f_pred_barQ.star_A <- function(QY.init, samp_data) {
        h_bars <- pred.hbars(samp_data, fit_h_reg_obj, NetInd_k)
        h_iptw <- h_bars$df_h_bar_vals$h
        determ.Q <- samp_data$determ.Q
        if (!is.na(coef(m.Q.star_h_A))) {
          off <- qlogis(QY.init)
          QY.star <- plogis(off + coef(m.Q.star_h_A)*h_iptw)
          #************************************************
          # print("determ.Q"); print(sum(determ.Q))
          #************************************************
          # QY.star[determ.Q] <- 1
          QY.star[determ.Q] <- samp_data[determ.Q, Ynode] # will be incorrect when W's are resampled
          #************************************************
          return(QY.star)
        } else {
          return(QY.init)  
        }
      }
      # TMLE B - Predict E[Y|g_star] (QY.star) for each i, based on the intercept epsilon update model for E[Y|C^Y] (m.Q.star_h_B)
      .f_pred_barQ.star_B <- function(QY.init, samp_data) {
        determ.Q <- samp_data$determ.Q
        if (!is.na(coef(m.Q.star_h_B))) {
          off <- qlogis(QY.init)
          QY.star_B <- plogis(off + coef(m.Q.star_h_B))
          #************************************************
          # QY.star[determ.Q] <- 1
          QY.star_B[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect when W's are resampled
          # print("QY.star_B"); print(QY.star_B)
          #************************************************
          return(QY.star_B)
        } else {
          return(QY.init)
        }
      }
      # iid TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the intercept epsilon update model for E[Y|C^Y] (m.Q.star_h_B)
      .f_pred_bariidQ.star_B <- function(iidQY.init, samp_data) {
        determ.Q <- samp_data$determ.Q
        if (!is.na(coef(m.iidQ.star_reg_B))) {
          iidoff <- qlogis(iidQY.init)
          iidQY.star_B <- plogis(iidoff + coef(m.iidQ.star_reg_B))
          #************************************************
          # QY.star[determ.Q] <- 1
          iidQY.star_B[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect when W's are resampled
          # print("iidQY.star_B"); print(iidQY.star_B)
          #************************************************
          return(iidQY.star_B)
        } else {
          return(iidQY.init)
        }
      }
  		# get an estimate of fi_W (hold ALL W's fixed at once) - a component of TMLE Var
  		.f.gen.fi_W <- function(NetInd_k, emp_netW) {
  			determ_df <-  data.frame(determ.g=data$determ.g, determ.Q=data$determ.Q)
        resamp_A <- f.gen.A.star(k, data.frame(emp_netW,subset(data, select=nFnode)), f.g.star, f.g_args)
  			samp_dataA <- .f.allCovars(k, NetInd_k, resamp_A, Anode)
  			resamp_A_fixW <- data.frame(emp_netW, samp_dataA, subset(data, select=c(nFnode, Ynode)), determ_df)
        # *******fi_W based on Q,N.init model ******
        QY.init_fixW <- .f_pred_barQ.init(m.Q.init, resamp_A_fixW)
        fi_W_init <- QY.init_fixW   # vers 2
        # *******fi_W based on Q,N.star models (A & B) ******
        if (alphaTMLE_B_only) {
          QY.star_fixW_A <- rep_len(0, nrow(resamp_A_fixW))
          QY.star_fixW_B <- .f_pred_barQ.star_B(QY.init_fixW, resamp_A_fixW)
        } else {
          QY.star_fixW_A <- .f_pred_barQ.star_A(QY.init_fixW, resamp_A_fixW)
          QY.star_fixW_B <- .f_pred_barQ.star_B(QY.init_fixW, resamp_A_fixW)          
        }
  			return(list(fi_W_init=fi_W_init, fi_W_star_A=QY.star_fixW_A, fi_W_star_B=QY.star_fixW_B))
  		}
  		# IPTW NETWORK TMLE
  		.f.gen_TMLEnetIPTW <- function(QY.init, samp_data, NetInd_k) {
  			g_iptw <- iptw_est(k=k, data=samp_data, node_l=node_l, m.gN=m.g0N, f.g.star=f.g.star, f.g_args=f.g_args, 
  								          family=family, NetInd_k=NetInd_k, lbound=lbound, max_npwt=max_npwt, f.g0=f.g0, f.g0_args=f.g0_args)
  			determ.Q <- samp_data$determ.Q
  			if (!is.na(coef(m.Q.star_iptw))) {
  				off <- qlogis(QY.init)
  				QY.star <- plogis(off + coef(m.Q.star_iptw)*g_iptw)
          #************************************************
  				# QY.star[determ.Q] <- 1
          QY.star[determ.Q] <- samp_data[determ.Q, Ynode]          
          #************************************************          
  				return(QY.star)
  			}
  		}
  		#-------------------------------------------		  			
  		# Main body of .f.gen.reps()
  		#-------------------------------------------
  		resamp_d <- .f.gen.sample(NetInd_k, n) # Get a random sample of all A and W
  		QY_gstar_mle <- .f_pred_barQ.init(m.Q.init, resamp_d) # QY.init (G-Comp estimator) - est probY based on model for Q_Y
      iidQY_gstar_mle <- .f_pred_barQ.init(m.iidQ.init, resamp_d) # QY.init (G-Comp estimator) - est probY based on model for Q_Y
  		#-------------------------------------------
      if (alphaTMLE_B_only) {
        QY_gstar_TMLE_A <- rep_len(0, n) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
        QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_iidTMLE_B <- .f_pred_bariidQ.star_B(iidQY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_TMLE_IPTW <- rep_len(0, n) # IPTW NETWORK TMLE
      } else {
        QY_gstar_TMLE_A <- .f_pred_barQ.star_A(QY_gstar_mle, resamp_d) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
        QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_iidTMLE_B <- .f_pred_bariidQ.star_B(iidQY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_TMLE_IPTW <- .f.gen_TMLEnetIPTW(QY_gstar_mle, resamp_d, NetInd_k) # IPTW NETWORK TMLE
      }
  		fi_Ws_list <- .f.gen.fi_W(NetInd_k, emp_netW) # Get fi_W - hold W fixed to observed values

      # print("nFriends under g_star"); print(table(resamp_d$nFriend)/ nrow(resamp_d))
      # print("Gcomp distr under g_star (MC)"); print(table(round(QY_gstar_mle,4)) / nrow(resamp_d))
      # print("MC gcomp"); print(length(QY_gstar_mle)); print(head(cbind(resamp_d,QY_gstar_mle), 200))
      # print("mean MC gcomp"); print(mean(QY_gstar_mle))

  		# Put all estimators together and add names (defined in G_D_W_1_nms outside of this function):
      mean_psis_all <- c(mean(QY_gstar_mle), mean(QY_gstar_TMLE_A), mean(QY_gstar_TMLE_B), mean(QY_gstar_iidTMLE_B), mean(QY_gstar_TMLE_IPTW), 
                        fi_Ws_list$fi_W_init, fi_Ws_list$fi_W_star_A, fi_Ws_list$fi_W_star_B)
  		names(mean_psis_all) <- G_D_W_1_nms
  		return(mean_psis_all)
  	} # end of .f.gen.reps()

  	#-------------------------------------------		  			
  	# Main body of .f_ests_all()
  	#-------------------------------------------
  	# Creating one data matrix to estimate all D*_W_i, for i=1,..,n
    # n_count <- n  	
  	# n_subsamp_vec <- NULL
  	# n_fullsamp_vec <- NULL
  	# idx_count <- 0
  	# idx_vec <- 1
  	# NetReInd_all <- NULL
  	# nF_vec <- NULL
  	# j_i_in_Fj_vec <- NULL
  	# j_fullsamp_vec <- NULL
  	# j_i_in_Fj_list <- list()
    # n_tot <- n
  	#-------------------------------------------
  	all_ests_reps <- t(sapply(seq(n_MCsims), .f.gen.reps, NetInd_k, emp_netW))
                        # idx_vec, 
                        # fixW=data[seq(n), node_l$Wnodes],
                        # NetReInd_all, nsubsamp=n_tot, 
                        # nF_vec, j_i_in_Fj_vec, j_fullsamp_vec))
  	return(all_ests_reps)
  }

  #---------------------------------------------------------------------------------
  # Main body of a fcn get.MCS_ests(): MC evalution of the estimators
  #---------------------------------------------------------------------------------
  # Names of all the estimators calculated during MC simulation:
  G_D_W_1_nms <- c("gcomp_mle", "tmle_A","tmle_B", "iid.tmle_B", "tmle_iptw",
                  paste("fWi_init_", c(1:n), sep = ""), 
                  paste("fWi_star_A_", c(1:n), sep = ""),
                  paste("fWi_star_B_", c(1:n), sep = ""))
  #  # returns list of indices (friends) that will be affected when W_i for i is held fixed
  # .f_get_arr_D_Wi <- function(i, NetInd_k) {
  # 	# Find all j's s.t. i \in F_j - evaluate for such j's, the rest of j's stay constant
  # 	j_i_in_Fj <- which(NetInd_k == i, arr.ind=TRUE)
  # 	j_i_in_Fj <- unique(c(i, as.vector(j_i_in_Fj[, 1])))
  # 	return(j_i_in_Fj)
  # }
  # # returns an array (or list) of indices for each i that needs to be resampled when W_i is fixed
  # # this includes not only all js s.t. i is in Fj but also entire network of each j
  # .f_get_arr_fullsamp <- function(i, NetInd_k, NetVec) {
  # 	j_i_in_Fj <- .f_get_arr_D_Wi(i, NetInd_k)
  # 	# Get j's & friends of j for final sample
  # 	j_samp <- unique(as.vector(c(j_i_in_Fj, as.vector(unlist(NetVec[j_i_in_Fj])))))
  # 	return(j_samp)
  # }
  # # list of all i's in F_j, by i=1,..,n
  # NetVec_k_D_Wi <- lapply(seq(n), .f_get_arr_D_Wi, NetInd_k) 
  # # list of full sample
  # NetVec_D_fullsamp <- lapply(seq(n), .f_get_arr_fullsamp, NetInd_k, NetVec) 
  #---------------------------------------------------------------------------------

  # Creating matrix of W's (fixed at observed Ws, for evalution of fi_W)
  netW <- NULL
  for (Wnode in node_l$Wnodes) {
    netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[,Wnode], Wnode)))
  }
  emp_netW <- netW

  # Allow this part to loop, until desired MCS prob_epsilon for all estimators is reached:
  nrepeat <- 1
  psis_reps <- NULL
  G_comp_D_star_W_reps <- NULL
  repeat {		  						  								
    G_comp_D_star_W_reps <- rbind(G_comp_D_star_W_reps, .f_ests_all(NetInd_k, emp_netW))
    # G_comp_D_star_W_reps <- rbind(G_comp_D_star_W_reps, .f_ests_all(NetInd_k, NetVec_k_D_Wi, NetVec_D_fullsamp))
  	psi_est_mean <- apply(G_comp_D_star_W_reps, 2, mean, na.rm = T)
  	psi_est_var <- apply(G_comp_D_star_W_reps, 2, var, na.rm = T)
  	psi_percerr <- 2 * abs(psi_est_mean * max.err_eps) # estimate the maximum allowed epsilon for each estimator, based pre-defined % error:  
  	# prob_epsilon <- psi_est_var / ((n_MCsims*nrepeat) * (max.err_eps)^2)
  	prob_percerr <- psi_est_var / ((n_MCsims*nrepeat) * (psi_percerr)^2)
  	prob_percerr[psi_est_var < 0.0001] <- 0.0001
  	fin_ests_sel <- c(1:3) # final vec of estimators for which error is measured
  	if ( (all(prob_percerr[fin_ests_sel] < 0.05)) | (nrepeat >= 100)) {
  		break
  	}
  	nrepeat <- nrepeat + 1
  }		
  # print("nrepeat"); print(nrepeat)
  return(psi_est_mean)
}

#---------------------------------------------------------------------------------
# USE THIS AS A TEMPLATE FOR REWRITE OF get.MCS_ests()
# G-Comp & TMLEs: Use Monte-Carlo to estimate psi under stochastic \bar{h}^* mixture density
#---------------------------------------------------------------------------------
  # For given data, take Q[Y|cY]=m.Q.init and calcualte est. of psi under \bar{h}^* using Monte-Carlo integration:
  # * Draw from the distributions of W and \bar{g}^*
  # * Recalculate Y^c under \bar{g}^*;
  # * Repeat nrep times and average.
get.MCS_ests_hstar <- function(data, data_net, MC_fit_params, fit_h_reg_obj, max.err_eps, family="binomial") {
  n <- nrow(data)
  k <- MC_fit_params$k
  # n_MCsims <- MC_fit_params$n_MCsims # not used

  lbound <- MC_fit_params$lbound
  node_l <- MC_fit_params$node_l
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
  Ynode <- node_l$Ynode
  iidW_flag <- MC_fit_params$iidW_flag
  NetInd_k <- MC_fit_params$NetInd_k
  # NetVec <- MC_fit_params$NetVec # not used

  m.Q.init <- MC_fit_params$m.Q.init
  m.Q.star_h_A <- MC_fit_params$m.Q.star_h_A
  m.Q.star_h_B <- MC_fit_params$m.Q.star_h_B
  f.g.star <- MC_fit_params$f.g.star; f.g_args <- MC_fit_params$f.g_args
  f.g0 <- MC_fit_params$f.g0; f.g0_args <- MC_fit_params$f.g0_args
  hstar <- MC_fit_params$hstar; hgN <- MC_fit_params$hgN; h_tilde <- MC_fit_params$h_tilde # not used

  # 02/07/12 Eliminated the need to make two passes through the data for monte-carlo # Creates one matrix of all Q(Y|W_i)
  .f.gen.reps <- function(nrep, NetInd_k, emp_netW, emp_netA)   {
    .f.gen.sample <- function(NetInd_k, n, useEmpD = TRUE) { # Generate full sample   c=(A,W) under g* (and iid or not for Ws)
      if (!useEmpD) {
        # ----------------------------
        # NEED TO SAMPLE FROM THE distribution of \bar{h}^* - need to sample c's under mixture distribution
        # -
        # <HERE>
        # # -
        # if (iidW_flag) {  # Version 1: Resample W's with replacement
        #   resamp_idx <- sample(c(1:n), n, replace=TRUE)
        #   netW <- NULL # get all network W's from the original data under resampled IDs
        #   for (Wnode in node_l$Wnodes) {
        #     netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[resamp_idx,Wnode], Wnode)))
        #   }
        #   full_dfW <- netW
        # }
        # resamp_A <- f.gen.A.star(k, full_dfW, f.g.star, f.g_args) # - this is wrong, as it gives a sample from g^*, not \bar{h}^c
        # full_dfA <- .f.allCovars(k, NetInd_k, resamp_A, Anode)
        # # ----------------------------
      } else {
        resamp_idx <- c(1:n)
        # full_dfW <- emp_netW
        # full_dfA <- emp_netA
      }
      determ_df <-  data.frame(determ.g=data$determ.g[resamp_idx], determ.Q=data$determ.Q[resamp_idx]) # get deterministic nodes, also resampled, since fcn of W's
      # Y_resamp <- subset(data, select=Ynode) # use observed Y's
      # nFriend <- subset(data, select=nFnode) # nFriends (network) is never resampled
      # resamp_data <- data.frame(full_dfW, full_dfA, Y_resamp, nFriend, determ_df)
      resamp_data <- data.frame(data_net, determ_df)
      return(resamp_data)
    }
    .f_pred_barQ.init <- function(samp_data) { # MLE - Predict E[Y|g_star] (QY.init) for each i, based on the initial model for E[Y|C^Y] (m.Q.init)
      determ.Q <- samp_data$determ.Q #deterministic nodes for Q
      QY <- predict(m.Q.init, newdata=samp_data, type="response")
      QY[determ.Q] <- samp_data[determ.Q, Ynode]
      return(QY)
    }
    .f_pred_barQ.star_A <- function(QY.init, samp_data) { # TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the coefficient epsilon update model for E[Y|C^Y] (m.Q.star_h_A)      
      h_bars <- pred.hbars(samp_data, fit_h_reg_obj, NetInd_k)
      h_iptw <- h_bars$df_h_bar_vals$h
      determ.Q <- samp_data$determ.Q
      if (!is.na(coef(m.Q.star_h_A))) {
        off <- qlogis(QY.init)
        QY.star <- plogis(off + coef(m.Q.star_h_A)*h_iptw)
        QY.star[determ.Q] <- samp_data[determ.Q, Ynode]
        return(QY.star)
      } else {
        return(QY.init)  
      }
    }
    .f_pred_barQ.star_B <- function(QY.init, samp_data) { # TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the intercept epsilon update model for E[Y|C^Y] (m.Q.star_h_B)
      determ.Q <- samp_data$determ.Q
      if (!is.na(coef(m.Q.star_h_B))) {
        off <- qlogis(QY.init)
        QY.star_B <- plogis(off + coef(m.Q.star_h_B))
        QY.star_B[determ.Q] <- samp_data[determ.Q, Ynode]
        return(QY.star_B)
      } else {
        return(QY.init)
      }
    }
    #-------------------------------------------
    resamp_d <- .f.gen.sample(NetInd_k=NetInd_k, n=n, useEmpD=TRUE) # Get a random sample of all A and W
    QY_gstar_mle <- .f_pred_barQ.init(resamp_d) # QY.init (G-Comp estimator) - est probY based on model for Q_Y
    QY_gstar_TMLE_A <- .f_pred_barQ.star_A(QY_gstar_mle, resamp_d) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
    QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
    #-------------------------------------------
    res <- c(gcomp_mle=sum(QY_gstar_mle*(h_tilde)), tmle_A=sum(QY_gstar_TMLE_A*(h_tilde)), tmle_B=sum(QY_gstar_TMLE_B*(h_tilde))) / sum((h_tilde))    
    # print("weighted results"); print(res)
    res
  }
  psi_est_mean <- .f.gen.reps(1, NetInd_k)
  return(psi_est_mean)
}

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
pred.hbars <- function(new_data=NULL, fit_h_reg_obj, NetInd_k) {
   pred_h_fcn <- function(logit_sep_k, m.gAi_vec) {
      .f_pred_h_k <- function(k, sel_k_indx, m.gAi_vec) {
        k_sel <- sel_k_indx
        A_nms_arr <- colnames(indA)[c(1:(k+1))]
        indA <- indA[,c(1:(k+1)), drop=FALSE]

        if (is.null(h_form)) {
          W_nms_arr <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
          W_nms_arr <- W_nms_arr[!is.na(W_nms_arr)]         
        } else {
          W_nms_arr <- all.vars(as.formula(h_form))[-1]
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
      if (logit_sep_k) {
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
    logit_sep_k <- fit_h_reg_obj$logit_sep_k
    node_l <- fit_h_reg_obj$node_l
    gform <- fit_h_reg_obj$gform
    h_form <- fit_h_reg_obj$h_form
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
    if (is.null(h_form)) {
      W_nms <- unlist(netW_namesl)
    } else {
      W_nms <- all.vars(as.formula(h_form))[-1]
    }
    if (!(node_l$nFnode%in%W_nms)) { W_nms <- c(W_nms, node_l$nFnode) }
    cY_mtx <- cbind(determ_cols, as.matrix(new_data[, c(node_l$nFnode, W_nms, netA_names)]))
    #---------------------------------------------------------------------
    # MAIN BODY OF THE FUNCTION
    if (h_user==FALSE) {
      P.hbar.c <- pred_h_fcn(logit_sep_k, fit_h_reg_obj$m.gAi_vec_g)
      P.hbar.star.c <- pred_h_fcn(logit_sep_k, fit_h_reg_obj$m.gAi_vec_gstar)
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

# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars <- function(data, h_fit_params) {
    #---------------------------------------------------------------------------------
    # 1) return observed network data cY_mtx if is.null(f.g_name)
    # 2) pass only one netW, which can be separate for g_star and g_0
    # SAMPLE A LARGE DATASET of cY's for given functions g* or g0, of size p*N for some p
    # Get rid of the loop, by assigning .f.allCovars ALL AT ONCE
    # NOTE (OS 06/02/15): The only reason we need to pass netW and netW_full is because g_0 and g^* 
    # are assumed to be based on the same sW! This shouldn't be the case, need to allow sW_g ad sW_gstar to be different
    #---------------------------------------------------------------------------------
    # get_hfit_data <- function(cY_mtx, k, Anode, NetInd_k, netW, f.g_name, f.g_args, p, misval = 0L)  {
    #   samp_g_data <- function(df_sel) {
    #     resamp_A <- f.gen.A.star(k, df_sel, f.g_name, f.g_args)
    #     resamp_netA <- .f.allCovars(k, NetInd_k, resamp_A, Anode, misval = misval)
    #     fit.g_data <- cbind(netW, resamp_netA)
    #     return(fit.g_data)
    #   }
    #   if (is.null(f.g_name)) { return(cY_mtx) }
    #   n <- nrow(netW)
    #   df_samp_g <- samp_g_data(netW)  # big resampled matrix of c's (repeated data p times)
    #   fit.g_data_large <- matrix(nrow=(n*p), ncol=ncol(df_samp_g))
    #   colnames(fit.g_data_large) <- colnames(df_samp_g)
    #   fit.g_data_large[c(1:n),] <- as.matrix(df_samp_g)
    #   for (i in (2:p)) {
    #     fit.g_data_large[c(((i - 1) * n + 1):(n * i)), ] <- as.matrix(samp_g_data(netW))
    #   }
    #   return(fit.g_data_large)
    # }
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
    fit_h_fcn <- function(fit.g_data, p, logit_sep_k, netW_namesl, indA) {
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
      if (is.null(h_form)) {
        stop("NOT IMPLEMENTED, SUPPLY h_form")
        # W_nms_arr <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
        # W_nms_arr <- W_nms_arr[!is.na(W_nms_arr)]
      } else {
        W_nms_arr <- all.vars(as.formula(h_form))[-1]
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
    k <- h_fit_params$k
    node_l <- h_fit_params$node_l
    NetInd_k <- h_fit_params$NetInd_k
    lbound <- h_fit_params$lbound
    max_npwt <- h_fit_params$max_npwt

    logit_sep_k=h_fit_params$logit_sep_k
    h_user=h_fit_params$h_user; h_user_fcn=h_fit_params$h_user_fcn; h_form=h_fit_params$h_form
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
    hform_covars <- all.vars(as.formula(h_form))[-1]
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
    if (is.null(h_form)) {
      W_nms <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
      W_nms <- W_nms[!is.na(W_nms)]
    } else {
      W_nms <- all.vars(as.formula(h_form))[-1]
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
    A_nms <- netvar2(Anode, c(0:k))
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
    P.hbar.c <- fit_h_fcn(fit.g_data = fit.g0_dat, p = p_h0, logit_sep_k = logit_sep_k, netW_namesl = netW_namesl, indA = indA)
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
    P.hbar.star.c <- fit_h_fcn(fit.g_data = fit.gstar_dat, p = n_samp_g0gstar, logit_sep_k = logit_sep_k, netW_namesl = netW_namesl, indA = indA)
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
                          k=k,
                          lbound=lbound,
                          logit_sep_k=logit_sep_k,
                          h_user=h_user,
                          h_user_fcn=h_user_fcn,
                          fit_fastmtx=fit_fastmtx,
                          node_l=node_l,
                          netW_namesl=netW_namesl,
                          gform = gform,
                          h_form = h_form,
                          netA_names=colnames(indA),
                          determ_cols_Friend=determ_cols_Friend,
                          determ_cols_fitted=determ_cols, 
                          cY_mtx_fitted=cY_mtx)
    return(list(df_h_bar_vals=df_h_bar_vals, fit_h_reg_obj=fit_h_reg_obj))
}

#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
#---------------------------------------------------------------------------------
.f_get_all_ests <- function(data, data_net, est_obj, MCeval_hstar) {
  n <- nrow(data)
  node_l <- est_obj$node_l
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
	Ynode <- node_l$Ynode
	Y <- data[, Ynode]
	determ.Q <- data[, "determ.Q"]
  determ.g <- data[, "determ.g"]

  #************************************************
  # 2 IID ESTIMATORS (completely ignore the network information)
  iidQY.init <- data[, "iidQY.init"]
  iidoff <- qlogis(iidQY.init)
  #************************************************
  # IID IPTW weights and estimator (ignores the network):
  #************************************************
  iidg <- iptw_est(k=est_obj$k, data=data, node_l=node_l, m.gN=est_obj$m.iidg0N, f.g.star=est_obj$f.g.star, f.g_args=est_obj$f.g_args, family=est_obj$family, 
                  NetInd_k=est_obj$NetInd_k, lbound=est_obj$lbound, max_npwt=est_obj$max_npwt, f.g0=est_obj$f.g0, f.g0_args=est_obj$f.g0_args, 
                  iidIPTW=TRUE)
  Y_iidIPTW <- Y
  Y_iidIPTW[!determ.Q] <- Y[!determ.Q] * iidg[!determ.Q]
  print("iid IPW Est:"); print(mean(Y_iidIPTW))
  #************************************************
  # IID TMLE based on weighted univariate update (espsilon is the intercept), IID IPTW weights and IID Q_N:
  #************************************************
  ctrl <- glm.control(trace=FALSE, maxit=1000)
  SuppressGivenWarnings(m.iidQ.star_reg_B <- glm(Y ~ offset(iidoff), data=data.frame(Y=data[, Ynode], iidoff=iidoff), weights=iidg,
                                            subset=!determ.Q, family=est_obj$family, control=ctrl), GetWarningsToSuppress(TRUE))
  iidQY.star_B <- Y
  if (!is.na(coef(m.iidQ.star_reg_B))) iidQY.star_B  <- plogis(iidoff + coef(m.iidQ.star_reg_B))

  #************************************************
  # NETWORK ESTIMATORS:
  QY.init <- data[, "QY.init"] # initial Q fit
  off <- qlogis(QY.init)  # offset

  #************************************************
  # IPTW_h estimator (based on h^*/h_N clever covariate):
  #************************************************
  fit.hbars_t <- system.time(h_bars <- fit.hbars(data=data, h_fit_params=est_obj)) # fit the clever covariat
  fit.hbars_t.new <- system.time(h_bars.new <- fit.hbars.new(data=data, h_fit_params=est_obj)) # fit the clever covariat
  df_h_bar_vals <- h_bars$df_h_bar_vals
  fit_h_reg_obj <- h_bars$fit_h_reg_obj
  h_iptw <- df_h_bar_vals$h
  Y_IPTW_h <- Y
  Y_IPTW_h[!determ.Q] <- Y[!determ.Q] * h_iptw[!determ.Q]
  # print("IPW Est (h)"); print(mean(Y_IPTW_h))

  print("time to fit h_bars"); print(fit.hbars_t)
  print("h est"); print(head(df_h_bar_vals))
  print("time to fit h_bars new"); print(fit.hbars_t.new)
  print("h est new"); print(head(h_bars.new$df_h_bar_vals))

  #************************************************
  # TMLE A: estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET
  #************************************************
  ctrl <- glm.control(trace=FALSE, maxit=1000)
  SuppressGivenWarnings(m.Q.star_reg_A <- glm(Y ~ -1 + h_iptw + offset(off), data=data.frame(Y=data[, Ynode], off=off, h_iptw=h_iptw),
                                						subset=!determ.Q, family=est_obj$family, control=ctrl), GetWarningsToSuppress(TRUE))
	QY.star <- Y
	if (!is.na(coef(m.Q.star_reg_A))) QY.star  <- plogis(off + coef(m.Q.star_reg_A) * h_iptw)
  #************************************************
  # TMLE B: estimate the TMLE update via weighted univariate ML (espsilon is intercept)
  #************************************************
  ctrl <- glm.control(trace=FALSE, maxit=1000)
  SuppressGivenWarnings(m.Q.star_reg_B <- glm(Y ~ offset(off), data=data.frame(Y=data[, Ynode], off=off), weights=h_iptw,
                                            subset=!determ.Q, family=est_obj$family, control=ctrl), GetWarningsToSuppress(TRUE))
  QY.star_B <- Y
  if (!is.na(coef(m.Q.star_reg_B))) QY.star_B  <- plogis(off + coef(m.Q.star_reg_B))
  #************************************************
  # IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
	# 02/16/13: IPTW estimator (Y_i * prod_{j \in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
	g_iptw <- iptw_est(k=est_obj$k, data=data, node_l=node_l, m.gN=est_obj$m.g0N, f.g.star=est_obj$f.g.star, f.g_args=est_obj$f.g_args, family=est_obj$family, 
                      NetInd_k=est_obj$NetInd_k, lbound=est_obj$lbound, max_npwt=est_obj$max_npwt, f.g0=est_obj$f.g0, f.g0_args=est_obj$f.g0_args)
  Y_IPTW_net <- Y
  Y_IPTW_net[!determ.Q] <- Y[!determ.Q] * g_iptw[!determ.Q]
  #************************************************
  # IPTW-based clever covariate TMLE (based on FULL likelihood factorization), covariate based fluctuation
  #************************************************
	SuppressGivenWarnings(m.Q.star_iptw <- glm(Y ~ -1 + g_iptw + offset(off),
                                						data=data.frame(Y=data[, Ynode], off=off, g_iptw=g_iptw),
                                						subset=!determ.Q, family=est_obj$family, control=ctrl),
                                						GetWarningsToSuppress(TRUE))

  parsubmodel_fits <- rbind(coef(m.Q.star_reg_A), coef(m.Q.star_reg_B), coef(m.Q.star_iptw))
  rownames(parsubmodel_fits) <- c("epsilon (covariate)", "alpha (intercept)", "iptw epsilon (covariate)")
  print("parsubmodel_fits"); print(parsubmodel_fits)

  #************************************************
  # Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
  # TO DO: create a MC specific object with defined structure: models on Q's, data info, etc...
  # TO DO: create a MC evaluation that can work for any g^* or \bar{g}^*...
	#************************************************
  # MC_sim_params <- list(k=k, node_l=node_l, NetInd_k=NetInd_k, lbound=lbound, max_npwt=max_npwt, n_MCsims=n_MCsims, alphaTMLE_B_only=alphaTMLE_B_only,
  #                       m.iidQ.init=m.iidQ.init, m.g0N=m.g0N, m.Q.init=m.Q.init,
  #                       m.Q.star_h_A=m.Q.star_reg_A, m.Q.star_h_B=m.Q.star_reg_B, m.Q.star_iptw=m.Q.star_iptw,
  #                       f.g.star=f.g.star, f.g_args=f.g_args, f.g0=f.g0, f.g0_args=f.g0_args,
  #                       hstar=df_h_bar_vals$h.star_c, hgN=df_h_bar_vals$h_c, h_tilde=df_h_bar_vals$h,
  #                       m.iidQ.star_reg_B=m.iidQ.star_reg_B)
  MC_fit_params <- append(est_obj,
                      list(m.iidQ.star_reg_B=m.iidQ.star_reg_B,
                          m.Q.star_h_A=m.Q.star_reg_A, 
                          m.Q.star_h_B=m.Q.star_reg_B, 
                          m.Q.star_iptw=m.Q.star_iptw, 
                          hstar=df_h_bar_vals$h.star_c,
                          hgN=df_h_bar_vals$h_c,
                          h_tilde=df_h_bar_vals$h))
  # run M.C. evaluation estimating psi under g^*:
  syst1 <- system.time(MCS_res <- get.MCS_ests(data=data, data_net=data_net,  MC_fit_params=MC_fit_params, fit_h_reg_obj=fit_h_reg_obj))
	# print("time to run MCS: "); print(syst1);
  #************************************************
  # TO DO: come up with a better way to handle various estimators below:
  psi_iid.tmle_B <- MCS_res[names(MCS_res)=="iid.tmle_B"]
  psi_mle <- MCS_res[names(MCS_res)=="gcomp_mle"]
  psi_tmle_A <- MCS_res[names(MCS_res)=="tmle_A"]
  psi_tmle_B <- MCS_res[names(MCS_res)=="tmle_B"]
  psi_tmle_iptw <- MCS_res[names(MCS_res)=="tmle_iptw"]
  psi_iptw_h <- mean(Y_IPTW_h)  # IPTW estimator based on h - clever covariate
  psi_iptw <- mean(Y_IPTW_net)  # IPTW estimator based on full g factorization (prod(g))
  psi_iid.iptw <- mean(Y_iidIPTW)  # IPTW estimator based on full g factorization (prod(g))
  # MCS_estimators <- MCS_res[names(MCS_res) %in% c("gcomp_mle", "tmle", "tmle_iptw")]

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

  # run the MC estimate evaluation by sampling from hstar instead of g^*
  # TO DO: remove this and use one M.C. sampling function that can take either g^* or bar{g}^*
  if (MCeval_hstar) { 
    syst2 <- system.time(MCS_res_hstar <- get.MCS_ests_hstar(data=data, data_net=data_net,  MC_fit_params=MC_sim_params,  fit_h_reg_obj=fit_h_reg_obj))
    noMC_psi_tmle_A <- MCS_res_hstar[names(MCS_res_hstar)=="tmle_A"]
    noMC_psi_tmle_B <- MCS_res_hstar[names(MCS_res_hstar)=="tmle_B"]
    noMC_psi_mle <- MCS_res_hstar[names(MCS_res_hstar)=="gcomp_mle"]
  } else {
    noMC_psi_tmle_A <- noMC_psi_tmle_B <- noMC_psi_mle <- 0
  }
	#-------------------------------------------
  return(list( iid.iptw = psi_iid.iptw,
               iid.tmle_B= psi_iid.tmle_B,
               tmle_A = psi_tmle_A,
               noMC_tmle_A = noMC_psi_tmle_A,
               tmle_B = psi_tmle_B,
               noMC_tmle_B = noMC_psi_tmle_B,
               tmle_iptw = psi_tmle_iptw,
               iptw_h = psi_iptw_h,
               iptw = psi_iptw,
               mle = psi_mle,
               noMC_mle = noMC_psi_mle,
               fWi_init_A=fWi_init_A, fWi_star_A=fWi_star_A,
               fWi_init_B=fWi_init_B, fWi_star_B=fWi_star_B,
               fWi_init_tmleiptw=fWi_init_tmleiptw,
               h_iptw=h_iptw, iptw_reg=g_iptw, QY.init=QY.init, QY.star=QY.star))
}

#---------------------------------------------------------------------------------
# MAIN TMLE ESTIMATOR FUNCTION
#---------------------------------------------------------------------------------
tmlenet <- function(data, Anode, AnodeDET = NULL, Wnodes, iidW_flag=FALSE, Ynode, YnodeDET=NULL, nFnode, 
                    k_max, IDnode=NULL, NETIDnode,
                    Qform=NULL, Q.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
                    gform=NULL, gbound=0.005, max_npwt=50, 
                    g.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
                    f.g1.star, f.g1_args, f.g2.star=NULL, f.g2_args=NULL, 
                    h_f.g0=NULL, h_f.g0_args=NULL, h_user_fcn = NULL, h_form = NULL,
                    h_logit_sep_k = FALSE, W_indep=FALSE,  family = "binomial", 
                    alpha  = 0.05, verbose=FALSE, n_MCsims=ceiling(sqrt(nrow(data))), n_samp_g0gstar=20, 
                    alphaTMLE_B_only = TRUE) {

  # if (is.na(h_form)) h_form <- NULL
  message("Running tmlenet with... "); 
  message("Qform: "%+%Qform); message("gform: "%+%gform); message("hform: "%+%h_form)

  # DEBUGGING:
  # memory & speed profiling
  # Rprof("tmle_run_memuse.out", memory.profiling=TRUE)

  #------------------------------------------------------------------------------
  # Which esimators to evaluate
  #------------------------------------------------------------------------------
  # alphaTMLE_B_only <- TRUE # if TRUE, only evalute the intercept TMLE (TMLE_B)
  MCeval_hstar <- FALSE   # if FALSE, do not evaluate MC wrt density h_star

  #------------------------------------------------------------------------------
  # MONTE-CARLO SIMULATION PARAMETERS
  #------------------------------------------------------------------------------
  n_MCsims <- as.integer(n_MCsims)  # number of times to sample for MC sim  
  max.err_est <- 0.1    # maximum percent error for MCS estimators

  #------------------------------------------------------------------------------
  # ESTIMATION OF h & h^* PARAMETERS
  #------------------------------------------------------------------------------
  n_samp_g0gstar <- n_samp_g0gstar # number of times to sample from g^* or g_0 (if known) for estimation of \bar{h}^* and \bar{h}_0
	logit_sep_k = h_logit_sep_k # Fit each k value separately (T) or all at once (F) (now passed as arguments to tmlenet)
  f.g0 <- h_f.g0
  f.g0_args <- h_f.g0_args
	d <- data
	n <- nrow(d)
	k <- k_max
  h_user <- !(is.null(h_user_fcn))

  # **********
  #Q: 02/15/14: IS IT A GOOD IDEA TO HAVE THE SAME (OR ANY) UPPER BOUND ON g_star? Doesn't make any sense...
  # **********
	if(length(gbound)==1) gbound <- c(gbound, 1-gbound)
 	#------------------------------------------------------------------------------
	# Netwk ids strings to list of friend indices (NetVec) and matrix of friend indices (NetInd_k)
  .f_getNetIndices <- function (d) {
    .f.mkvecNetID <- function(Net_str) lapply(Net_str, function(Net_str_i) unlist(strsplit(Net_str_i, ' ', fixed=TRUE)))
      # Netwk ids strings to arr of ID vectors (filled up with trailing NA's)
    .f.mkarrNetID <- function(Net_str) t(sapply(Net_str, function(Net_str_i) {
                      netwk <- unlist(strsplit(Net_str_i, ' ',
                      fixed=TRUE))
                      return(c(netwk, rep_len(NA,k-length(netwk))))
                      } ))
    NetIDVec <- .f.mkvecNetID(d[,NETIDnode])  # get list of network ID vectors from network strings for each i
    NetVec <- lapply(NetIDVec, function(NetID) as.numeric(sapply(NetID, function(x) which(x==as.vector(d[,IDnode]))))) # convert into list of network row #s
    NetID_k <- .f.mkarrNetID(d[,NETIDnode]) # get array of network IDs
    NetInd_k <- apply(NetID_k, 2, function(k_ID) {  # make NetID_k into array of network indices (row #s)
                  sapply(as.vector(k_ID), function(x) {
                    if (is.na(x)) {
                      NA
                    } else { 
                      which(x==as.vector(d[,IDnode]))
                    }
                   })
                })
    return(list(NetVec=NetVec, NetInd_k=NetInd_k))
  }
  NetInd_l <- .f_getNetIndices(d)
  NetVec <- NetInd_l$NetVec
  NetInd_k  <- NetInd_l$NetInd_k

  #------------------------------------------------------------------------------
  # List of node names and networks of W's and A's (netA and netW)
  #------------------------------------------------------------------------------
  node_l <- list(IDnode=IDnode, Anode=Anode, Wnodes=Wnodes, Ynode=Ynode, nFnode=nFnode, NETIDnode=NETIDnode)

  if (is.null(AnodeDET)) determ.g <- rep_len(FALSE,n) else 
      determ.g <- (d[, AnodeDET] == 1)
  if (is.null(YnodeDET)) determ.Q <- rep_len(FALSE,n) else 
      determ.Q <- (d[, YnodeDET] == 1)

	netW <- NULL
	for (Wnode in node_l$Wnodes) {
	 	netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, d[,Wnode], Wnode)))
	}
	netA <- data.frame(.f.allCovars(k, NetInd_k, d[,node_l$Anode], node_l$Anode))
	df_AllWs <- netW	 
	net_d <- cbind(ID=d[, node_l$IDnode], netW, netA, subset(d, select=c(node_l$nFnode,node_l$Ynode)))
	# print("net_d"); print(head(net_d, 10))

	#-------------------------------------------					
	# Fit initial model for Q(Y|A,W) under oberved data (m.Q.init) - move this to .f_get_all_ests() to avoid confusion
	#-------------------------------------------  
  print("fitting the inital Q"); print(Qform)
	m.Q.init <- .f.est(net_d[!determ.Q,], Qform, family=family)  # Set Y=Y_i when determ.Q_i==1
  QY.init <- d[, Ynode] # set deterministic nodes
  QY.init[!determ.Q] <- predict(m.Q.init, newdata=net_d[!determ.Q,], type="response") # predict p(Y) for non determ nodes 

  iidQform <- formula(x="Y~"%+%Anode%+%"+"%+%paste(Wnodes, collapse="+"))
  print("fitting the inital iid Q"); print(iidQform)
  m.iidQ.init <- .f.est(net_d[!determ.Q,], iidQform, family=family)  # Set Y=Y_i when determ.Q_i==1
  iidQY.init <- d[, Ynode] # set deterministic nodes
  iidQY.init[!determ.Q] <- predict(m.iidQ.init, newdata=net_d[!determ.Q,], type="response") # predict p(Y) for non determ nodes 

	d_sel <- data.frame(subset(d, select=unlist(node_l)), determ.g=determ.g, determ.Q=determ.Q, QY.init=QY.init, iidQY.init=iidQY.init)
  # node_l <- c(node_l, gform=gform, Qform=Qform, iidW_flag=iidW_flag)
  # reg_models <- list(gform=gform, Qform=Qform, iidW_flag=iidW_flag)
  # print("m.Q.init"); print(summary(m.Q.init)); print("QY.init fit"); print(QY.init);

	#-------------------------------------------							
	# Fit the model for g_0(A,W) - move this to .f_get_all_ests() to avoid confusion
	#-------------------------------------------
	m.g0N <- .f.est(net_d[!determ.g,], gform, family=family) # Set A=0 when determ.g==1

  iidgform <- Anode%+%"~"%+%paste(Wnodes, collapse="+")
  print("Using the formula to fit iid g: "); print(iidgform)
  m.iidg0N <- .f.est(net_d[!determ.g,], iidgform, family=family) # Set A=0 when determ.g==1

  #-------------------------------------------              
  # Create an object with model estimates, data & network information that is passed on to estimation procedure
  #-------------------------------------------
  # 1) define parameters for MC estimation of the substitution estimators
  # 2) define parameters for estimation of the efficient weights h(A^s|W^s)
  est_obj <- list(
    k=k, 
    node_l=node_l,
    NetInd_k=NetInd_k,
    lbound=gbound[1], 
    alphaTMLE_B_only=alphaTMLE_B_only, n_MCsims=n_MCsims, 
    max.err_eps=max.err_est,  # error tolerance for the mean/var M.C. estimatef
    max_npwt=max_npwt,  # NOT IMPLEMENTED YET: cap the prop weights scaled at max_npwt (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000)
    family = family,
    m.iidg0N=m.iidg0N, m.iidQ.init=m.iidQ.init,
    m.g0N=m.g0N, m.Q.init=m.Q.init,
    f.g0=f.g0, f.g0_args=f.g0_args,
    logit_sep_k=logit_sep_k, 
    h_user=h_user, 
    h_user_fcn=h_user_fcn, 
    Qform = Qform,
    gform = gform,
    h_form = h_form,
    iidW_flag = iidW_flag,
    n_samp_g0gstar=n_samp_g0gstar)
	
  #-------------------------------------------							
	# Run TMLE for each g.star and/or ATE
	#-------------------------------------------
  est_obj_g1 <- append(est_obj, list(f.g.star=f.g1.star, f.g_args=f.g1_args))
	tmle_g1_out <- .f_get_all_ests(data=d_sel, data_net=net_d, est_obj=est_obj_g1, MCeval_hstar=MCeval_hstar)
	if (!is.null(f.g2.star)) {
    est_obj_g2 <- append(est_obj, list(f.g.star=f.g2.star, f.g_args=f.g2_args))
		tmle_g2_out <- .f_get_all_ests(data=d_sel, data_net=net_d, est_obj=est_obj_g2, MCeval_hstar=MCeval_hstar)
	}
	else {
		tmle_g2_out <- NULL
	}
	#-------------------------------------------
	# Estimate IC-based Variance (as. Var) and CIs (based on f_W and fY)
	#-------------------------------------------
	# use estimates of fWi (hold all W's fixed at once),
	# loop over all intersecting friends networks
	# calculate R_ij*(fW_i*fW_j) - see page 33 Network paper vdL
	#----------------------------------------------------------------------------------
	# Helper function to calculate cross product sum of correlated f_Wi (see p.33 of vdL)
  # New fast method for as var calculation (matrix vs)
  .f_est_sigmas_new <- function(QY.init, QY.star, fWi_A, fWi_B, fWi_tmleiptw, h_iptw, est_iptw_h, iptw_reg, est_iptw, alphaTMLE_B_only) {
      var_tmle_A_iidIC <- var_tmle_B_iidIC <- var_iid.tmle_B_iidIC <- var_tmleiptw_iidIC_1stO <- var_tmleiptw_iidIC_2ndO <- var_iptw_h_iidIC <- var_iptw_iidIC_1stO <- var_iptw_iidIC_2ndO <- var_tmle_A_Q.init <- var_tmle_B_Q.init <- 0
      # get the connectivity n_by_n mtx (1 indicates intersection of friendship sets)
      # returns 1) 1st order and 2) 1st and 2nd order connections
      get.Fiintersectmtx <- function() {
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
      connectmtx_obj <- get.Fiintersectmtx()
      connectmtx_1stO <- connectmtx_obj$conn_ind_mtx_1st
      connectmtx_2ndO <- connectmtx_obj$conn_ind_mtx_2ndO

      # TMLE A (clever covariate update): Inference based on the iid IC analogy, QY.init - initial Q model predictions, h_iptw - h_tilde
      if (!alphaTMLE_B_only) {
        iidIC_tmle_A <- h_iptw*(d_sel[,node_l$Ynode]-QY.init) + fWi_A
        var_tmle_A_iidIC <- est.sigma_fsum(get.crossprodmtx(iidIC_tmle_A), connectmtx_1stO)
      }

      # **NEW** TMLE B (weights model update): Inference based on the iid IC analogy:
      iidIC_tmle_B <- h_iptw*(d_sel[,node_l$Ynode]-QY.init) + fWi_B
      var_tmle_B_iidIC <- est.sigma_fsum(get.crossprodmtx(iidIC_tmle_B), connectmtx_1stO)
      # simple iid estimator of the asymptotic variance (no adjustment for correlated observations i,j):
      var_iid.tmle_B_iidIC <- mean((iidIC_tmle_B)^2)

      # TMLE based on iptw clever covariate (more non-parametric)
      if (!alphaTMLE_B_only) {
        iidIC_tmleiptw <- iptw_reg*(d_sel[,node_l$Ynode]-QY.init) + fWi_tmleiptw
        var_tmleiptw_iidIC_1stO <- est.sigma_fsum(get.crossprodmtx(iidIC_tmleiptw), connectmtx_1stO)
        var_tmleiptw_iidIC_2ndO <- est.sigma_fsum(get.crossprodmtx(iidIC_tmleiptw), connectmtx_2ndO)
      }

      # **NEW** IPTW based on the mixture density clever covariate (h): 
      iidIC_iptw_h <- h_iptw*(d_sel[,node_l$Ynode]) - est_iptw_h
      var_iptw_h_iidIC <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw_h), connectmtx_1stO)

      # IPTW:
      iidIC_iptw <- iptw_reg*(d_sel[,node_l$Ynode]) - est_iptw
      var_iptw_iidIC_1stO <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw), connectmtx_1stO)
      var_iptw_iidIC_2ndO <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw), connectmtx_2ndO)
 
      # Inference based on the EIC, with factorization into orthogonal components sigma2_DY and sigma2_W_N
      # sigma2_DY_i are independent (since they are conditioned on W,A)
      # sigma2_W_N_i are dependent => need to take double sum of their crossprod among dependent units
      if (!alphaTMLE_B_only) {
        D_star_Yi.Qinit <- h_iptw * (d_sel[,node_l$Ynode] - QY.init) # h*(Y-Q_bar_N):
        sigma2_DY <- (1/n) * sum(D_star_Yi.Qinit^2)  # Sum_{i} (D_star_Yi)^2

        fW_A_crossprod <- get.crossprodmtx(fWi_A)
        sigma2_W_N_A <- est.sigma_fsum(fW_A_crossprod, connectmtx_1stO)
        var_tmle_A_Q.init <- sigma2_W_N_A + sigma2_DY

        # **NEW** TMLE B (weights model update)
        fW_B_crossprod <- get.crossprodmtx(fWi_B)
        sigma2_W_N_B <- est.sigma_fsum(fW_B_crossprod, connectmtx_1stO)
        var_tmle_B_Q.init <- sigma2_W_N_B + sigma2_DY

        # D_star_Yi.Qstar <- h_iptw * (d_sel[,node_l$Ynode] - QY.star)
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
      return(list(var_tmle_A_iidIC=abs(var_tmle_A_iidIC), # new
                  var_tmle_B_iidIC=abs(var_tmle_B_iidIC), # new
                  var_iid.tmle_B_iidIC=abs(var_iid.tmle_B_iidIC), # no adjustment for correlations i,j
                  var_tmleiptw_iidIC_1stO=abs(var_tmleiptw_iidIC_1stO),
                  var_tmleiptw_iidIC_2ndO=abs(var_tmleiptw_iidIC_2ndO),
                  var_iptw_h_iidIC=abs(var_iptw_h_iidIC), # new
                  var_iptw_iidIC_1stO=abs(var_iptw_iidIC_1stO),
                  var_iptw_iidIC_2ndO=abs(var_iptw_iidIC_2ndO),
                  var_tmle_A_Q.init=abs(var_tmle_A_Q.init),
                  var_tmle_B_Q.init=abs(var_tmle_B_Q.init)
                  # var_tmle_A_Q.star_cons=abs(var_tmle_A_Q.star_cons)
              ))
  }

  # create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
  .f_make_EYg_obj <- function(tmle_g_out, tmle_g2_out=NULL) {
    # get estimates of as. var and CIs:
    if (is.null(tmle_g2_out)) {
      tmle_g2_out <- list()
      tmle_g2_out$QY.star <- tmle_g2_out$fWi_init_A <- tmle_g2_out$fWi_init_B <- tmle_g2_out$fWi_star_A <- tmle_g2_out$fWi_star_B <- tmle_g2_out$fWi_init_tmleiptw <- tmle_g2_out$h_iptw <- tmle_g2_out$iptw_reg <- 0
      tmle_g2_out$tmle_A <- tmle_g2_out$tmle_B <- tmle_g2_out$iid.tmle_B <- tmle_g2_out$tmle_iptw <- tmle_g2_out$iptw_h <- tmle_g2_out$iptw <- tmle_g2_out$iid.iptw <- tmle_g2_out$mle <- 0
      tmle_g2_out$noMC_tmle_A <- tmle_g2_out$noMC_tmle_B <- tmle_g2_out$noMC_mle <- 0
    }

    tmle_A <- tmle_g_out$tmle_A - tmle_g2_out$tmle_A
    noMC_tmle_A <- tmle_g_out$noMC_tmle_A - tmle_g2_out$noMC_tmle_A

    tmle_B <- tmle_g_out$tmle_B - tmle_g2_out$tmle_B
    iid.tmle_B <- tmle_g_out$iid.tmle_B - tmle_g2_out$iid.tmle_B
    noMC_tmle_B <- tmle_g_out$noMC_tmle_B - tmle_g2_out$noMC_tmle_B

    tmle_iptw <- tmle_g_out$tmle_iptw - tmle_g2_out$tmle_iptw
    iptw_h <- tmle_g_out$iptw_h - tmle_g2_out$iptw_h
    iptw <- tmle_g_out$iptw - tmle_g2_out$iptw
    iid.iptw <- tmle_g_out$iid.iptw - tmle_g2_out$iid.iptw
    
    mle = tmle_g_out$mle - tmle_g2_out$mle
    noMC_mle = tmle_g_out$noMC_mle - tmle_g2_out$noMC_mle
 
    EY_g.star <- list(  tmle_A = as.vector(tmle_A),
                        noMC_tmle_A = as.vector(noMC_tmle_A),
                        tmle_B = as.vector(tmle_B),
                        iid.tmle_B = as.vector(iid.tmle_B),
                        noMC_tmle_B = as.vector(noMC_tmle_B),
                        tmle_iptw = as.vector(tmle_iptw),
                        iptw_h = as.vector(iptw_h),
                        iptw = as.vector(iptw),
                        iid.iptw = as.vector(iid.iptw),
                        mle = as.vector(mle),
                        noMC_mle = as.vector(noMC_mle)
                        )

    # IS THIS INCORRECT???? SHOULD THIS BE A DIFFERENCE FOR 2 g^*'s? OR THEY WILL BE THE SAME? 
    # THEY WILL BE THE SAME, AS THEY DON'T DEPEND on g^*...
    QY.init <- tmle_g_out$QY.init
    QY.star <- tmle_g_out$QY.star

    #***** which fWi should we use? Based on mQ.init or mQ.star? - mQ.star, because fWi has psi over i component?
    fWi_init_A <- tmle_g_out$fWi_init_A - tmle_g2_out$fWi_init_A  # tmle_A (update model with clever covar)
    fWi_init_B <- tmle_g_out$fWi_init_B - tmle_g2_out$fWi_init_B  # tmle_B (update model with weights)
    fWi_star_A <- tmle_g_out$fWi_star_A - tmle_g2_out$fWi_star_A
    fWi_star_B <- tmle_g_out$fWi_star_B - tmle_g2_out$fWi_star_B
    # fWi <- tmle_g_out$fWi_star - tmle_g2_out$fWi_star
    #*****
    fWi_init_tmleiptw <- tmle_g_out$fWi_init_tmleiptw - tmle_g2_out$fWi_init_tmleiptw  # tmle_iptw (update model with iptw clever covar)
    h_iptw <- tmle_g_out$h_iptw - tmle_g2_out$h_iptw
    iptw_reg <- tmle_g_out$iptw_reg - tmle_g2_out$iptw_reg

    ICnms <- c("tmle_A_iidIC","tmle_B_iidIC","iid.tmle_B_iidIC","tmleiptw_iidIC_1stO", "tmleiptw_iidIC_2ndO","iptw_h_iidIC", "iptw_iidIC_1stO", "iptw_iidIC_2ndO", "tmle_A_Q.init", "tmle_B_Q.init")
    varsnms <- "var_"%+%ICnms
    CIsnms <- "CI_"%+%ICnms

    .f_gen_var_CI <- function(est_IC_name, psi) {
      # get an est. of CI:
      .f_est_CI <- function(sigma2_N, psi_tmle) {
        z_alpha <- qnorm(1-alpha/2)
        CI_tmle <- c(psi_tmle - z_alpha*sqrt(sigma2_N)/sqrt(n), 
                      psi_tmle + z_alpha*sqrt(sigma2_N)/sqrt(n))
        names(CI_tmle) <- c("LBCI", "UBCI")
        return(CI_tmle)
      }
      var_nm <- paste0("var_", est_IC_name)
      CI_nm <- paste0("CI_", est_IC_name)
      res <- list(as.vector(tmle_vars_obj[[var_nm]]/n), .f_est_CI(tmle_vars_obj[[var_nm]], psi))
      names(res) <- c(var_nm, CI_nm)
      return(res)
    }

    # estiamte asymptotic variance for each estimator
    tmle_vars_obj <- .f_est_sigmas_new(QY.init=QY.init, QY.star=QY.star, fWi_A=fWi_init_A, fWi_B=fWi_init_B, 
                                        fWi_tmleiptw=fWi_init_tmleiptw, h_iptw=h_iptw, est_iptw_h=iptw_h, 
                                        iptw_reg=iptw_reg, est_iptw=iptw, alphaTMLE_B_only=alphaTMLE_B_only)

    EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmle_A_iidIC",tmle_A))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmle_B_iidIC",tmle_B)) # new (weight-based update model)
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("iid.tmle_B_iidIC",tmle_B))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmleiptw_iidIC_1stO",tmle_iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmleiptw_iidIC_2ndO",tmle_iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("iptw_h_iidIC",iptw_h)) # new (clever covariate-based iptw)
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("iptw_iidIC_1stO",iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("iptw_iidIC_2ndO",iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmle_A_Q.init",tmle_A))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmle_B_Q.init",tmle_B))
    # EY_g.star <- c(EY_g.star, .f_gen_var_CI("tmle_A_Q.star_cons",tmle_A))
    return(EY_g.star)
  }

  EY_g1.star <- .f_make_EYg_obj(tmle_g_out=tmle_g1_out)

	if (!is.null(f.g2.star)) {	
    EY_g2.star <- .f_make_EYg_obj(tmle_g_out=tmle_g2_out)
    ATE <- .f_make_EYg_obj(tmle_g_out=tmle_g1_out, tmle_g2_out=tmle_g2_out)
	} else {
		EY_g2.star <- NULL
		ATE <- NULL
	}
  # print("EY_g1.star"); print(EY_g1.star)
  # print("EY_g2.star"); print(EY_g2.star)
	estimates <- list(EY_g1.star=EY_g1.star, EY_g2.star=EY_g2.star, ATE=ATE)
	tmlenet <- list(estimates=estimates)
	class(tmlenet) <- "tmlenet"	

  out_sim <- rbind(EY_g1.star$tmle_A, EY_g1.star$tmle_B, EY_g1.star$iid.tmle_B, EY_g1.star$tmle_iptw, EY_g1.star$iptw_h, EY_g1.star$iptw, EY_g1.star$iid.iptw, EY_g1.star$mle)
  rownames(out_sim) <- c("tmle_A", "tmle_B", "iid.tmle_B", "tmle_iptw", "iptw_h", "iptw", "iid.iptw", "mle")
  print("Estimates w/ MC eval:"); print(out_sim)

  noMC_out_sim <- rbind(EY_g1.star$noMC_tmle_A, EY_g1.star$noMC_tmle_B, EY_g1.star$noMC_mle)
  rownames(noMC_out_sim) <- c("noMC_tmle_A","noMC_tmle_B","noMC_mle")
  # print("Estimates w/out MC eval:"); print(noMC_out_sim)

  # DEBUGGING:
  # Rprof(NULL)
	return(tmlenet)
}
#----------------------------------------------------------------------------------	  
# #   browser {base}
# R Documentation
# Environment Browser
# Description
# Interrupt the execution of an expression and allow the inspection of the environment where browser was called from.
#----------------------------------------------------------------------------------	  