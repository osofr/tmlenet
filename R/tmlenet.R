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
  #   Kmax - Constant for maximum number of friends any observation can have (cannot be exceeded)
  #   IDnode - Optional subject identifier variable in data, if not supplied the network in NETIDnode is assumed to be indexed by row numbers
  #   NETIDnode - String variable in data identifying subject's network by IDs or row numbers (space separated)
  #   Qform -  Regression formula for outcome, Q 
  #   gform -  Regression formula for treatment mechanism, g
  #   hform - Regression formula for the clever covariate, P(A^s | W^s) / P(A^{*s} | W^s) used for estimating under g_0 and g^*.
  #   AnodeDET - Column name for indicators of deterministic values of A, coded as (TRUE/FALSE) or (1/0) (to be removed: TRUE sets A=0), (to add: the observations with AnodeDET=TRUE/1 are assumed to have constant value for their Anode)
  # (TO ADD)  AnodeDETfun - function that evaluates to TRUE for observations with deterministically assigned A (alternative to AnodeDET)
  #   YnodeDET - Column name for indicators of deterministic values of Y, coded as (TRUE/FALSE) or (1/0) (to be removed: TRUE sets Y=1), (to add: the observations with YnodeDET=TRUE/1 are assumed to have constant value for their Ynode)
  # (TO ADD)  YnodeDETfun - function that evaluates to TRUE for observations with deterministically assigned Y (alternative to YnodeDET)
  #   Q.SL.library - SuperLearner libraries for outcome, Q (NOT IMPLEMENTED)
  #   g.SL.library - SuperLearner libraries for treatment mechanism, g (NOT IMPLEMENTED)
  #   gbound - One value for symmetrical bounds on g(A|W), or a vector containing upper and lower bounds
  #   f_g1star - Function for dynamic treatment regimen, can take any variables in data as arguments
  #   args_f_g1star - Additional arguments to be passed to f_g1star
  #   f_g2star - Optional dynamic treatment regimen for contrasts
  #   args_f_g2star - Optional additional arguments to be passed to f_g2star
  #   h_f.g0 - When known, a function for generating A under true treatment mechanism, g0. Used only during estimation of the clevel covariate, h, via logistic regression after sampling A from this function
  #   h_f.g0_args - Additional arguments to be passed to h_f.g0
  #   h_user - Should a user supplied function be used to calculate the clever covariate, h
  #   h_user_fcn - User supplied function to calculate the clever covariate, h
  #   h_logit_sep_k - Flag for fitting a separate logistic regression for each strata of nFnode, used during estimation of the clever covariate, h
  #   family - Family specification for regression models, defaults to binomial. CURRENTLY ONLY BINOMIAL FAMILY IS IMPLEMENTED.
  #   alpha - alpha-level for CI calculation
  #   verbose - Flag for controlling printing of messages (NOT IMPLEMENTED)
  #   nQ.MCsims=ceiling(sqrt(nrow(data))) - number of Monte-Carlo simulation to perform for evaluation of psi under g^*
  #   ng.MCsims=20 - number of Monte-Carlo simulations to perform for evaluation of \h^* under g^*
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
  # ATE - additive treatment effect for EY_g1.star-EY_g0 (default) or EY_g1.star-EY_g2.star (if f_g2star is provided)
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
# G-Comp & TMLEs: Use Monte-Carlo to estimate psi under stochastic g^* 
#---------------------------------------------------------------------------------
  # For given data, take Q[Y|cY]=m.Q.init and calcualte est. of psi under g*=f.g.star using Monte-Carlo integration:
  # * W_i can be iid or not (in latter case W are not resampled);
  # * Draw from the distributions of W and g*(A|W), keeping N_i constant, recalculate cY and cA each time;
  # * Recalculate Y^c under g^*
  # * Repeat nrep times and average.
#---------------------------------------------------------------------------------
get.MCS_ests <- function(data, MC_fit_params, fit_h_reg_obj) {
  family <- MC_fit_params$family

  alphaTMLE_B_only <- MC_fit_params$alphaTMLE_B_only
  n <- nrow(data)
  k <- MC_fit_params$Kmax
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
  f.g.star <- MC_fit_params$f.g.star
  f.g_args <- MC_fit_params$f.g_args
  f.g0 <- MC_fit_params$f.g0
  f.g0_args <- MC_fit_params$f.g0_args

  # 02/07/12 Eliminated the need to make two passes through the data for monte-carlo # Creates one matrix of all Q(Y|W_i)
  .f_ests_all <- function(NetInd_k, emp_netW) {
    .f.gen.reps <- function(nrep, NetInd_k, emp_netW) 	{
      # Generate full sample   c=(A,W) under g* (and iid or not for Ws)
       .f.gen.sample <- function(NetInd_k, n) {
        if (iidW_flag) {  # Version 1: Resample W's with replacement
          resamp_idx <- sample(c(1:n), n, replace=TRUE)
          netW <- NULL # get all network W's from the original data under resampled IDs
          for (Wnode in node_l$Wnodes) {
            netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[resamp_idx,Wnode], Wnode)))
          }
          full_dfW <- netW
        } else {
        # Version 2: No resampling of W's (W's are indep. but not iid, using NPMLE that puts mass 1 on all obs i=1,...,N)
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
        h_bars <- pred.hbars.old(samp_data, fit_h_reg_obj, NetInd_k)
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
  		#-------------------------------------------
      if (alphaTMLE_B_only) {
        QY_gstar_TMLE_A <- rep_len(0, n) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
        QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_TMLE_IPTW <- rep_len(0, n) # IPTW NETWORK TMLE
      } else {
        QY_gstar_TMLE_A <- .f_pred_barQ.star_A(QY_gstar_mle, resamp_d) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
        QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_TMLE_IPTW <- .f.gen_TMLEnetIPTW(QY_gstar_mle, resamp_d, NetInd_k) # IPTW NETWORK TMLE
      }
  		fi_Ws_list <- .f.gen.fi_W(NetInd_k, emp_netW) # Get fi_W - hold W fixed to observed values

  		# Put all estimators together and add names (defined in G_D_W_1_nms outside of this function):
      mean_psis_all <- c(mean(QY_gstar_mle), mean(QY_gstar_TMLE_A), mean(QY_gstar_TMLE_B), mean(QY_gstar_TMLE_IPTW), 
                        fi_Ws_list$fi_W_init, fi_Ws_list$fi_W_star_A, fi_Ws_list$fi_W_star_B)
  		names(mean_psis_all) <- G_D_W_1_nms
  		return(mean_psis_all)
  	} # end of .f.gen.reps()

  	#-------------------------------------------
  	# Main body of .f_ests_all()
  	#-------------------------------------------
  	all_ests_reps <- t(sapply(seq(n_MCsims), .f.gen.reps, NetInd_k, emp_netW))
  	return(all_ests_reps)
  }
  #---------------------------------------------------------------------------------
  # Main body of a fcn get.MCS_ests(): MC evalution of the estimators
  #---------------------------------------------------------------------------------
  # Names of all the estimators calculated during MC simulation:
  G_D_W_1_nms <- c("gcomp_mle", "tmle_A","tmle_B", "tmle_iptw",
                  paste("fWi_init_", c(1:n), sep = ""),
                  paste("fWi_star_A_", c(1:n), sep = ""),
                  paste("fWi_star_B_", c(1:n), sep = ""))
  
  #---------------------------------------------------------------------------------
  # Creating matrix of W's (fixed at observed Ws, for evalution of fi_W)
  netW <- NULL
  for (Wnode in node_l$Wnodes) {
    netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[,Wnode], Wnode)))
  }
  emp_netW <- netW

  #---------------------------------------------------------------------------------
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
get.MCS_ests_hstar <- function(data, MC_fit_params, fit_h_reg_obj, max.err_eps, family="binomial") {
  n <- nrow(data)
  k <- MC_fit_params$Kmax
  lbound <- MC_fit_params$lbound
  node_l <- MC_fit_params$node_l
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
  Ynode <- node_l$Ynode
  iidW_flag <- MC_fit_params$iidW_flag
  NetInd_k <- MC_fit_params$NetInd_k
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
      }
      determ_df <-  data.frame(determ.g=data$determ.g[resamp_idx], determ.Q=data$determ.Q[resamp_idx]) # get deterministic nodes, also resampled, since fcn of W's
      resamp_data <- data.frame(data_net, determ_df)
      return(resamp_data)
    }
    .f_pred_barQ.init <- function(samp_data) { # MLE - Predict E[Y|g_star] (QY.init) for each i, based on the initial model for E[Y|C^Y] (m.Q.init)
    }
    .f_pred_barQ.star_A <- function(QY.init, samp_data) { # TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the coefficient epsilon update model for E[Y|C^Y] (m.Q.star_h_A)      
    }
    .f_pred_barQ.star_B <- function(QY.init, samp_data) { # TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the intercept epsilon update model for E[Y|C^Y] (m.Q.star_h_B)
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
  # m.Q.init$reg                # regression class (Yreg)

  nodes <- datNetObs$nodes
  Y <- datNetObs$noNA.Ynodevals # actual observed Y's
  determ.Q <- datNetObs$det.Y
  m.Q.init <- est_params_list$m.Q.init

  QY.init <- datNetObs$noNA.Ynodevals # getting all node vals, inc. deterministic
  QY.init[!datNetObs$det.Y] <- m.Q.init$getprobA1[!datNetObs$det.Y] # getting predictions P(Y=1) for non-DET Y
  off <- qlogis(QY.init)  # offset

  #************************************************
  # ESTIMATORS
  #************************************************

  #************************************************
  # h^*/h_N clever covariate:
  #************************************************
  fit.hbars_t.new <- system.time(h_bars.new <- fit.hbars(datNetObs = datNetObs, est_params_list = est_params_list)) # fit the clever covariat
  df_h_bar_vals <- h_bars$df_h_bar_vals
  fit_h_reg_obj <- h_bars$fit_h_reg_obj
  h_iptw <- df_h_bar_vals$h
  print("time to fit h_bars"); print(fit.hbars_t)
  print("h est"); print(head(df_h_bar_vals))

  #************************************************
  # IPTW_h estimator:
  #************************************************  
  Y_IPTW_h <- Y
  Y_IPTW_h[!determ.Q] <- Y[!determ.Q] * h_iptw[!determ.Q]
  # print("IPW Est (h)"); print(mean(Y_IPTW_h))

  # print("time to fit h_bars new"); print(fit.hbars_t.new)
  # print("h est new"); print(head(h_bars.new$df_h_bar_vals))

  #************************************************
  # TMLE A: estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET
  # #todo 19 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star_reg_A
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.star_reg_A <- glm(Y ~ -1 + h_iptw + offset(off), data = data.frame(Y = Y, off = off, h_iptw = h_iptw),
                                						subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))
	QY.star <- Y
	if (!is.na(coef(m.Q.star_reg_A))) QY.star <- plogis(off + coef(m.Q.star_reg_A) * h_iptw)

  #************************************************
  # TMLE B: estimate the TMLE update via weighted univariate ML (espsilon is intercept)
  # #todo 20 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star_reg_B
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.star_reg_B <- glm(Y ~ offset(off), data = data.frame(Y = Y, off = off), weights = h_iptw,
                                            subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))
  QY.star_B <- Y
  if (!is.na(coef(m.Q.star_reg_B))) QY.star_B <- plogis(off + coef(m.Q.star_reg_B))

  #************************************************
  # IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
	# 02/16/13: IPTW estimator (Y_i * prod_{j \in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
	# g_iptw <- iptw_est(k = est_params_list$Kmax, data = data, node_l = nodes, m.gN = est_params_list$m.g0N,
  #                      f.g.star = est_params_list$f.g.star, f.g_args = est_params_list$f.g_args, family = "binomial",
  #                      NetInd_k = est_params_list$NetInd_k, lbound = est_params_list$lbound, max_npwt = est_params_list$max_npwt,
  #                      f.g0 = est_params_list$f.g0, f.g0_args = est_params_list$f.g0_args)
  #  Y_IPTW_net <- Y
  #  Y_IPTW_net[!determ.Q] <- Y[!determ.Q] * g_iptw[!determ.Q]
  Y_IPTW_net <- rep(0, length(Y))

  #************************************************
  # IPTW-based clever covariate TMLE (based on FULL likelihood factorization), covariate based fluctuation
  # #todo 21 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star_iptw
  #************************************************
	# SuppressGivenWarnings(m.Q.star_iptw <- glm(Y ~ -1 + g_iptw + offset(off),
  #                                						data = data.frame(Y = Y, off = off, g_iptw = g_iptw),
  #                                						subset = !determ.Q, family = "quasibinomial", control = ctrl),
  #                                						GetWarningsToSuppress(TRUE))

  parsubmodel_fits <- rbind(coef(m.Q.star_reg_A), coef(m.Q.star_reg_B), coef(m.Q.star_iptw))
  rownames(parsubmodel_fits) <- c("epsilon (covariate)", "alpha (intercept)", "iptw epsilon (covariate)")
  print("parsubmodel_fits: "); print(parsubmodel_fits)




  #************************************************
  # Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
  # TO DO: create a MC specific object with defined structure: models on Q's, data info, etc...
  # TO DO: create a MC evaluation that can work for any g^*
	#************************************************
  MC_fit_params <- append(est_params_list,
                      list(m.Q.star_h_A = m.Q.star_reg_A,
                          m.Q.star_h_B = m.Q.star_reg_B,
                          m.Q.star_iptw = m.Q.star_iptw,
                          hstar = df_h_bar_vals$h.star_c,
                          hgN = df_h_bar_vals$h_c,
                          h_tilde = df_h_bar_vals$h))

  # run M.C. evaluation estimating psi under g^*:
  syst1 <- system.time(MCS_res <- get.MCS_ests(datNetObs = datNetObs,  MC_fit_params = MC_fit_params, fit_h_reg_obj = fit_h_reg_obj))
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

# MAIN TMLE ESTIMATOR FUNCTION:
# todo 24 (tmlenet) +0: need to sort out how to store sW.gstar and sW.g0 and how to pass those to fit.hbars
tmlenet <- function(data, Kmax, Anode, AnodeDET = NULL, Wnodes, Ynode, YnodeDET = NULL,
                    nFnode = "nF", IDnode = NULL, NETIDnode,
                    # CAN RETIRE Qform, gform & hform -> fully replaced by summary measures:
                    Qform = NULL, gform = NULL, hform = NULL,
                    f_g1star, args_f_g1star, f_g2star = NULL, args_f_g2star = NULL,
                    family = "binomial", alphaTMLE_B_only = TRUE,
                    gbound = 0.005, nQ.MCsims = ceiling(sqrt(nrow(data))), ng.MCsims = 20,
                    verbose = FALSE,
                    # can be either character or class %in% DefineEval.sVar
                    # can do S3 dispatch, when character just calls constructors for DefineEval.sVar
                    sW.g0, sW.gstar = NULL, sA) {

  #----------------------------------------------------------------------------------
  # ADDITIONAL ARGUMENTS (Removed from input args of tmlenet())
  #----------------------------------------------------------------------------------
  # alphaTMLE_B_only <- TRUE # if TRUE, only evalute the intercept TMLE (TMLE_B)
  iidW_flag <- FALSE
  Q.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction")
  g.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction")
  max_npwt <- 50
  h_logit_sep_k <- FALSE
  alpha <- 0.05
  h_f.g0 <- NULL
  h_f.g0_args <- NULL
  h_user_fcn <- NULL

  #----------------------------------------------------------------------------------
  # MONTE-CARLO SIMULATION PARAMETERS
  #----------------------------------------------------------------------------------
  nQ.MCsims <- as.integer(nQ.MCsims)  # number of times to sample for MC sim
  max.err_est <- 0.1    # maximum percent error for MCS estimators

  #----------------------------------------------------------------------------------
  # Input checks:
  # #todo 8 (tmlenet) +0: check all names exist in data (Anode, Wnodes, Ynode, etc...)
  # #todo 10 (tmlenet) +0: need default setting for sW.g0, sW.gstar & sA when these are left unscpecified
  # e.g., use dat.netW for sW and dat.netA for sA
  # #todo 11 (tmlenet) +0: Wnodes & Anode are no longer needed when sW, sA, give a warning that Wnodes,Anodes will be ignored
  #----------------------------------------------------------------------------------
  assert_that(is.data.frame(data))
  assert_that(is.count(Kmax))
  assert_that(is.count(nQ.MCsims))
  assert_that(is.count(ng.MCsims))
  # assert_that(is.integerish(Kmax))
  Kmax <- as.integer(Kmax)
  # d <- data
  # k <- Kmax
  nobs <- nrow(data)
  node_l <- list(IDnode = IDnode, Anode = Anode, Wnodes = Wnodes, Ynode = Ynode, nFnode = nFnode, NETIDnode = NETIDnode)

  # Defining Deterministic Y and A node flags:
  if (is.null(AnodeDET)) {determ.g <- rep_len(FALSE, nobs)} else {determ.g <- (data[, AnodeDET] == 1)}
  if (is.null(YnodeDET)) {determ.Q <- rep_len(FALSE, nobs)} else {determ.Q <- (data[, YnodeDET] == 1)}

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
  # Parse and evaluate the summary measures (in class DefineEval.sVar):
  # #todo 16 (tmletnet) +0: make below part of DatNet evaluation
  # #todo 17 (tmletnet) +5: pre-evaluate summary measures on small batch of data to get dims of sA & sW and to check for errors
  # Would need to sort out how to define NetInd_k for a subset of full data?
  # (NOT DONE) Check all sVar names are unique
  # (NOT DONE) If no sVars were defined (default), use netW and netA for sVars
  # (NOT DONE) Also accept sVars as character vectors / list, not only DefineEval.sVar
  #----------------------------------------------------------------------------------
  if (is.null(sW.gstar)) {
    sW.gstar <- sW.g0$clone
    sW.gstar$type <- "sW.gstar"
  }

  # MAKING THIS A TEST FOR EVALUATION OF sum. measures.
  # MOVING THE ACTUAL parsing steps to datnetW$make.sVar and datnetA$make.sVar
  testm.sW.g0 <- sW.g0$get.mat.sVar(data.df = data, netind_cl = netind_cl, addnFnode = node_l$nFnode)
  print("testm.sW.g0"); print(head(testm.sW.g0))
  testm.sW.gstar <- sW.gstar$get.mat.sVar(data.df = data, netind_cl  = netind_cl, addnFnode = node_l$nFnode)
  print("testm.sW.gstar"); print(head(testm.sW.gstar))
  testm.sA <- sA$get.mat.sVar(data.df = data, netind_cl = netind_cl)
  print("testm.sA"); print(head(testm.sA))

  # When sW / sA are just lists of character vectors need to capture the calling env and call constructors manually:
    # user.env <- parent.frame()
    # user.env_l <- list(user.env = user.env)
    # sW.g0 <- do.call(DefineEval.sVar$new, c(sW.g0, list(type = "sW.g0"), user.env_l))
    # sW.gstar <- do.call(DefineEval.sVar$new, c(sW.gstar, list(type = "sW.gstar"), user.env_l))
    # sA <- do.call(DefineEval.sVar$new, c(sA, list(type = "sA"), user.env_l))
  
  #----------------------------------------------------------------------------------
  # PARAMETERS FOR ESTIMATING h_g0 & h_gstar
  #----------------------------------------------------------------------------------
  f.g0 <- h_f.g0
  f.g0_args <- h_f.g0_args
  h_user <- !(is.null(h_user_fcn))

  #Q: 02/15/14: IS IT A GOOD IDEA TO HAVE THE SAME (OR ANY) UPPER BOUND ON g_star?
	if(length(gbound)==1) gbound <- c(gbound, 1 - gbound)

  #---------------------------------------------------------------------------------
  # BUILDING OBSERVED sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
  #---------------------------------------------------------------------------------
  datnetW <- DatNet$new(netind_cl = netind_cl, nodes = node_l, VarNodes = node_l$Wnodes, addnFnode = TRUE, misValRepl = TRUE)
  datnetW$make.sVar(Odata = data, sVar.object = sW.g0) # obsdat.sW <- datnetW$dat.sVar # datnetW$names.sVar

  print("head(obsdat.sW) as dat.sVar"); print(head(datnetW$dat.sVar))
  datnetW$fixmiss_sVar()
  print("head(obsdat.sW) after fixmiss_sVar:"); print(head(datnetW$dat.sVar))

  # EVERYTHING FAILS WHEN NOT REPLACING NA netVars WITH 0, SINCE MOST sA END UP AS NA:
  datnetA <- DatNet$new(netind_cl = netind_cl, nodes = node_l, VarNodes = node_l$Anode, misValRepl = TRUE)
  # datnetA <- DatNet$new(netind_cl = netind_cl, nodes = node_l, VarNodes = node_l$Anode)  
  datnetA$make.sVar(Odata = data, sVar.object = sA)
  print("head(obsdat.sA) as dat.sVar"); print(head(datnetA$dat.sVar))
  # datnetA$fixmiss_sVar()
  # print("head(obsdat.sA)after fixmiss_sVar: "); print(head(datnetA$dat.sVar))

  # IV) (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sA:
  # message("cancelled adding DET nodes to sVar since all sVar are automatically get added to A ~ predictors...")
  # obsdat.sW <- O.datnetW$add_deterministic(Odata = data, userDETcol = "determ.g")$dat.sVar
  # print(head(obsdat.sW))


  Yvals <- data[,node_l$Ynode]
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = Yvals, det.Y = determ.Q)$make.dat.sWsA()
  print("head(datNetObs$dat.sWsA)"); print(head(datNetObs$dat.sWsA))

  # browser()
  # stop()
  # Testing NA for visible det.Y and true observed Y as protected:
  # determ.Q <- c(FALSE, FALSE, FALSE, rep(TRUE, length(determ.Q)-3))
  # length(determ.Q) == length(Yvals)
  # datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = Yvals, det.Y = determ.Q)$make.dat.sWsA()

  #-----------------------------------------------------------
  # Defining and fitting regression for Y ~ sW + sA:
  #-----------------------------------------------------------
  # subset_expr <- try(parse(text="!misfun("%+%node_l$Ynode%+%")")[[1]])
  Yreg <- RegressionClass$new(outvar = node_l$Ynode,
                              predvars = c(datnetA$names.sVar, datnetW$names.sVar),
                              # predvars = c("netW3_sum", "sum_1mAW2_nets"),
                              subset = !determ.Q)
  m.Q.init <- BinOutModel$new(glm = FALSE, reg = Yreg)$fit(data = datNetObs)$predict()

  # datNetObs$YnodeVals       # visible Y's with NA for det.Y
  # datNetObs$det.Y           # TRUE/FALSE for deterministic Y's
  # datNetObs$noNA.Ynodevals  # actual observed Y's

  # m.Q.init$getoutvarnm      # reg outvar name (Ynode)
  # m.Q.init$getoutvarval     # original Yvals after setting det.Y obs to NA
  # m.Q.init$getprobA1        # predictions (for non-DET Y)
  # m.Q.init$getsubset        # valid subset (!det.Y)
  # m.Q.init$reg              # regression class (Yreg)

  # NO LONGER NEEDED HERE / REMOVE:
  QY.init <- datNetObs$YnodeVals # setting deterministic nodes
  QY.init[!datNetObs$det.Y] <- m.Q.init$getprobA1[!datNetObs$det.Y] # getting predictions P(Y=1)

  # print("summary(m.Q.init): "); print(summary(m.Q.init))
  print("coef(m.Q.init): "); print(coef(m.Q.init))

	#----------------------------------------------------------------------------------
  # Defining and fitting regression for A ~ sW:
	# Fit the model for g_0(A,W) - move this to get_all_ests() to avoid confusion
	#----------------------------------------------------------------------------------
  Areg <- RegressionClass$new(outvar = node_l$Anode,
                              predvars = c(datnetW$names.sVar),
                              # predvars = c("W1", "netW1_sum", "netW2_sum", "netW3_sum", "nFriends"),
                              subset = !determ.g)

  m.g0N <- BinOutModel$new(glm = FALSE, reg = Areg)$fit(data = datNetObs)

  # print("summary(m.g0N): "); print(summary(m.g0N))
  print("coef(m.g0N): "); print(coef(m.g0N))

  #----------------------------------------------------------------------------------
  # Create an object with model estimates, data & network information that is passed on to estimation procedure
  #----------------------------------------------------------------------------------
  # 1) define parameters for MC estimation of the substitution estimators
  # 2) define parameters for estimation of the efficient weights h(A^s|W^s)
  est_obj <- list(
                  alphaTMLE_B_only = alphaTMLE_B_only,
                  Kmax = Kmax, # REMOVE, already saved in DatNet
                  node_l = node_l, #REMOVE, already saved in DatNet
                  NetInd_k = NetInd_k, #REMOVE, already saved in DatNet

                  lbound = gbound[1],
                  max.err_eps = max.err_est,  # error tolerance for the mean/var M.C. estimate

                  # NOT IMPLEMENTED YET:
                  # cap the prop weights scaled at max_npwt
                  # (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000)
                  max_npwt = max_npwt,

                  family = family, #REMOVE, already saved in DatNet regs
                  m.g0N = m.g0N, m.Q.init = m.Q.init, # changed to BinOutModel objects

                  f.g0 = f.g0, f.g0_args = f.g0_args,
                  h_logit_sep_k = h_logit_sep_k,
                  h_user = h_user,
                  h_user_fcn = h_user_fcn,
                  gform = gform, #REMOVE, no longer used
                  hform = hform, #REMOVE, no longer used
                  iidW_flag = iidW_flag, #REMOVE?

                  n_MCsims = nQ.MCsims,        # REMOVE, keeping for now for compatibility
                  n_samp_g0gstar = ng.MCsims,  # REMOVE, keeping for now for compatibility
                  nQ.MCsims = nQ.MCsims,
                  ng.MCsims = ng.MCsims
                  )

  est_obj_g1 <- append(est_obj,
                      list(
                        f.g.star = f_g1star,
                        f.g_args = args_f_g1star)
                      )
  if (!is.null(f_g2star)) {
    est_obj_g2 <- append(est_obj,
                      list(
                        f.g.star = f_g2star,
                        f.g_args = args_f_g2star)
                      )
  }


  #----------------------------------------------------------------------------------
  # DEPRECATED... All regs are now defined via sW, sA summary measures...
  # Create net_d for fitting m.Q.init, m.g0N and m.h_g0, m.h_gstar
  # Check all RHS var names in Qform, gform, hform exist (NOT DONE)
  #----------------------------------------------------------------------------------
  if (is.null(Qform)) {
    # default to main terms in sW.g0 and sA
    Qform <- node_l$Ynode %+% " ~ " %+%
              paste0(colnames(datNetObs$datnetW$names.sVar), collapse="+") %+%
              "+" %+%
              paste0(colnames(datNetObs$datnetA$names.sVar), collapse="+")
  }
  if (is.null(gform)) {
    gform <- node_l$Anode %+% " ~ " %+% paste0(colnames(datNetObs$datnetW$names.sVar), collapse="+") # default to main terms in datNetObs$datnetW
  }
  if (is.null(hform)) {
    hform <- "s.A ~ " %+% paste0(colnames(datNetObs$datnetW$names.sVar), collapse="+") # default to main terms in datNetObs$datnetW
  }
  message("Running tmlenet with... ");
  message("Qform: " %+% Qform)
  message("gform: " %+% gform)
  message("hform: " %+% hform)

  net_d <- cbind(datNetObs$dat.sWsA, subset(data, select = node_l$Ynode))
  print("head(net_d)"); print(head(net_d, 5))
  #----------------------------------------------------------------------------------
  # (DEPRECATED) Defining and fitting regression for Y ~ sW + sA: Define df d_sel
  #----------------------------------------------------------------------------------
  m.Q.init.old <- f_est(net_d[!determ.Q,], Qform, family = family)
  QY.init.old <- data[, node_l$Ynode] # setting deterministic nodes
  QY.init.old[!determ.Q] <- predict(m.Q.init.old, newdata = net_d[!determ.Q,], type = "response") # predict p(Y) for non determ nodes
  print("coef(m.Q.init.old): "); print(coef(m.Q.init.old))
  m.g0N.old <- f_est(net_d[!determ.g,], gform, family = family) # Set A = 0 when determ.g == 1
  print("coef(m.g0N.old)"); print(coef(m.g0N.old))
  d_sel <- data.frame(subset(data, select = unlist(node_l)), determ.g = determ.g, determ.Q = determ.Q, QY.init = QY.init.old) # (DEPRECATED, TO BE REMOVED)
  print("head(d_sel) old: "); print(head(d_sel))

  #----------------------------------------------------------------------------------
	# DEPRECATED: Run TMLE univariate fluctuations for each g.star and/or ATE:
	#----------------------------------------------------------------------------------
  est_obj_g1$m.g0N <- m.g0N.old
  est_obj_g1$m.Q.init <- m.Q.init.old
  print(class(est_obj_g1$m.g0N))
  print(class(est_obj_g1$m.Q.init))

	tmle_g1_out <- get_all_ests.old(data = d_sel, est_obj = est_obj_g1)

  if (!is.null(f_g2star)) {
    est_obj_g2$m.g0N <- m.g0N.old
    est_obj_g2$m.Q.init <- m.Q.init.old

    tmle_g2_out <- get_all_ests.old(data = d_sel, est_obj = est_obj_g2)
	} else {
		tmle_g2_out <- NULL
	}

  #----------------------------------------------------------------------------------
  # ****** NEW ******
  # Run MC evaluation for substitution TMLE ests
  #----------------------------------------------------------------------------------
    est_obj_g1$m.g0N <- m.g0N
    est_obj_g1$m.Q.init <- m.Q.init

    tmle_g1_out <- get_all_ests(datNetObs = datNetObs, est_params_list = est_obj_g1)

    if (!is.null(f_g2star)) {
      est_obj_g2$m.g0N <- m.g0N
      est_obj_g2$m.Q.init <- m.Q.init

      tmle_g2_out <- get_all_ests(datNetObs = datNetObs, est_params_list = est_obj_g2)
    } else {
      tmle_g2_out <- NULL
    }

	#----------------------------------------------------------------------------------
	# Estimate IC-based Variance (as. Var) and CIs (based on f_W and fY)
	#----------------------------------------------------------------------------------
	# use estimates of fWi (hold all W's fixed at once),
	# loop over all intersecting friends networks
	# calculate R_ij*(fW_i*fW_j) - see page 33 Network paper vdL
	#----------------------------------------------------------------------------------
	# Helper function to calculate cross product sum of correlated f_Wi (see p.33 of vdL)
  # New fast method for as var calculation (matrix vs)
  .f_est_sigmas_new <- function(n, QY.init, QY.star, fWi_A, fWi_B, fWi_tmleiptw, h_iptw, est_iptw_h, iptw_reg, est_iptw, alphaTMLE_B_only) {
      var_tmle_A_iidIC <- var_tmle_B_iidIC <- var_iid.tmle_B_iidIC <- var_tmleiptw_iidIC_1stO <- var_tmleiptw_iidIC_2ndO <- var_iptw_h_iidIC <- var_iptw_iidIC_1stO <- var_iptw_iidIC_2ndO <- var_tmle_A_Q.init <- var_tmle_B_Q.init <- 0
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

      # TMLE A (clever covariate update): Inference based on the iid IC analogy, QY.init - initial Q model predictions, h_iptw - h_tilde
      if (!alphaTMLE_B_only) {
        iidIC_tmle_A <- h_iptw * (d_sel[,node_l$Ynode]-QY.init) + fWi_A
        var_tmle_A_iidIC <- est.sigma_fsum(get.crossprodmtx(iidIC_tmle_A), connectmtx_1stO)
      }

      # **NEW** TMLE B (weights model update): Inference based on the iid IC analogy:
      iidIC_tmle_B <- h_iptw * (d_sel[,node_l$Ynode]-QY.init) + fWi_B
      var_tmle_B_iidIC <- est.sigma_fsum(get.crossprodmtx(iidIC_tmle_B), connectmtx_1stO)
      # simple iid estimator of the asymptotic variance (no adjustment for correlated observations i,j):
      var_iid.tmle_B_iidIC <- mean((iidIC_tmle_B)^2)

      # TMLE based on iptw clever covariate (more non-parametric)
      if (!alphaTMLE_B_only) {
        iidIC_tmleiptw <- iptw_reg * (d_sel[,node_l$Ynode]-QY.init) + fWi_tmleiptw
        var_tmleiptw_iidIC_1stO <- est.sigma_fsum(get.crossprodmtx(iidIC_tmleiptw), connectmtx_1stO)
        var_tmleiptw_iidIC_2ndO <- est.sigma_fsum(get.crossprodmtx(iidIC_tmleiptw), connectmtx_2ndO)
      }

      # **NEW** IPTW based on the mixture density clever covariate (h): 
      iidIC_iptw_h <- h_iptw * (d_sel[,node_l$Ynode]) - est_iptw_h
      var_iptw_h_iidIC <- est.sigma_fsum(get.crossprodmtx(iidIC_iptw_h), connectmtx_1stO)

      # IPTW:
      iidIC_iptw <- iptw_reg * (d_sel[,node_l$Ynode]) - est_iptw
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
      return(list(
                  var_tmle_A_iidIC=abs(var_tmle_A_iidIC), # new
                  var_tmle_B_iidIC=abs(var_tmle_B_iidIC), # new
                  var_iid.tmle_B_iidIC=abs(var_iid.tmle_B_iidIC), # no adjustment for correlations i,j
                  var_tmleiptw_iidIC_1stO=abs(var_tmleiptw_iidIC_1stO),
                  var_tmleiptw_iidIC_2ndO=abs(var_tmleiptw_iidIC_2ndO),
                  var_iptw_h_iidIC=abs(var_iptw_h_iidIC), # new
                  var_iptw_iidIC_1stO=abs(var_iptw_iidIC_1stO),
                  var_iptw_iidIC_2ndO=abs(var_iptw_iidIC_2ndO),
                  var_tmle_A_Q.init=abs(var_tmle_A_Q.init),
                  var_tmle_B_Q.init=abs(var_tmle_B_Q.init)
              ))
  }

  # create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
  .f_make_EYg_obj <- function(tmle_g_out, tmle_g2_out=NULL) {
    # get estimates of as. var and CIs:
    if (is.null(tmle_g2_out)) {
      tmle_g2_out <- list()
      tmle_g2_out$QY.star <- tmle_g2_out$fWi_init_A <- tmle_g2_out$fWi_init_B <- tmle_g2_out$fWi_star_A <- tmle_g2_out$fWi_star_B <- tmle_g2_out$fWi_init_tmleiptw <- tmle_g2_out$h_iptw <- tmle_g2_out$iptw_reg <- 0
      tmle_g2_out$tmle_A <- tmle_g2_out$tmle_B <- tmle_g2_out$tmle_iptw <- tmle_g2_out$iptw_h <- tmle_g2_out$iptw <- tmle_g2_out$mle <- 0
    }

    tmle_A <- tmle_g_out$tmle_A - tmle_g2_out$tmle_A
    tmle_B <- tmle_g_out$tmle_B - tmle_g2_out$tmle_B
    tmle_iptw <- tmle_g_out$tmle_iptw - tmle_g2_out$tmle_iptw
    iptw_h <- tmle_g_out$iptw_h - tmle_g2_out$iptw_h
    iptw <- tmle_g_out$iptw - tmle_g2_out$iptw
    mle <- tmle_g_out$mle - tmle_g2_out$mle

    EY_g.star <- list(  tmle_A = as.vector(tmle_A),
                        tmle_B = as.vector(tmle_B),
                        tmle_iptw = as.vector(tmle_iptw),
                        iptw_h = as.vector(iptw_h),
                        iptw = as.vector(iptw),
                        mle = as.vector(mle)
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

    # tmle_iptw (update model with iptw clever covar):
    fWi_init_tmleiptw <- tmle_g_out$fWi_init_tmleiptw - tmle_g2_out$fWi_init_tmleiptw
    h_iptw <- tmle_g_out$h_iptw - tmle_g2_out$h_iptw
    iptw_reg <- tmle_g_out$iptw_reg - tmle_g2_out$iptw_reg

    ICnms <- c("tmle_A_iidIC",
              "tmle_B_iidIC",
              "iid.tmle_B_iidIC",
              "tmleiptw_iidIC_1stO",
              "tmleiptw_iidIC_2ndO",
              "iptw_h_iidIC",
              "iptw_iidIC_1stO",
              "iptw_iidIC_2ndO",
              "tmle_A_Q.init",
              "tmle_B_Q.init")
    # varsnms <- "var_" %+% ICnms
    # CIsnms <- "CI_" %+% ICnms

    .f_gen_var_CI <- function(n, est_IC_name, psi) {
      # get an est. of CI:
      .f_est_CI <- function(sigma2_N, psi_tmle) {
        z_alpha <- qnorm(1-alpha/2)
        CI_tmle <- c(psi_tmle - z_alpha*sqrt(sigma2_N) / sqrt(n),
                      psi_tmle + z_alpha*sqrt(sigma2_N) / sqrt(n))
        names(CI_tmle) <- c("LBCI", "UBCI")
        return(CI_tmle)
      }
      var_nm <- paste0("var_", est_IC_name)
      CI_nm <- paste0("CI_", est_IC_name)
      res <- list(as.vector(tmle_vars_obj[[var_nm]] / n), .f_est_CI(tmle_vars_obj[[var_nm]], psi))
      names(res) <- c(var_nm, CI_nm)
      return(res)
    }

    # estiamte asymptotic variance for each estimator
    tmle_vars_obj <- .f_est_sigmas_new(n = nobs, QY.init=QY.init, QY.star=QY.star, fWi_A=fWi_init_A, fWi_B=fWi_init_B,
                                        fWi_tmleiptw=fWi_init_tmleiptw, h_iptw=h_iptw, est_iptw_h=iptw_h,
                                        iptw_reg=iptw_reg, est_iptw=iptw, alphaTMLE_B_only=alphaTMLE_B_only)

    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "tmle_A_iidIC", tmle_A))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "tmle_B_iidIC", tmle_B)) # new (weight-based update model)
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "iid.tmle_B_iidIC", tmle_B))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "tmleiptw_iidIC_1stO", tmle_iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "tmleiptw_iidIC_2ndO", tmle_iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "iptw_h_iidIC", iptw_h)) # new (clever covariate-based iptw)
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "iptw_iidIC_1stO", iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "iptw_iidIC_2ndO", iptw))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "tmle_A_Q.init", tmle_A))
    EY_g.star <- c(EY_g.star, .f_gen_var_CI(n = nobs, "tmle_B_Q.init", tmle_B))
    return(EY_g.star)
  }

  EY_g1.star <- .f_make_EYg_obj(tmle_g_out=tmle_g1_out)

	if (!is.null(f_g2star)) {
    EY_g2.star <- .f_make_EYg_obj(tmle_g_out=tmle_g2_out)
    ATE <- .f_make_EYg_obj(tmle_g_out=tmle_g1_out, tmle_g2_out=tmle_g2_out)
	} else {
		EY_g2.star <- NULL
		ATE <- NULL
	}
	estimates <- list(EY_g1.star = EY_g1.star, EY_g2.star = EY_g2.star, ATE = ATE)
	tmlenet <- list(estimates = estimates)
	class(tmlenet) <- "tmlenet"

  out_sim <- rbind(EY_g1.star$tmle_A, EY_g1.star$tmle_B, EY_g1.star$tmle_iptw, EY_g1.star$iptw_h, EY_g1.star$iptw, EY_g1.star$mle)
  rownames(out_sim) <- c("tmle_A", "tmle_B", "tmle_iptw", "iptw_h", "iptw", "mle")
  print("Estimates w/ MC eval:"); print(out_sim)

	return(tmlenet)
}