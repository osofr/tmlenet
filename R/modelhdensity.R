#-----------------------------------------------------------------------------
# Fit and Predict the IPTW (clever covariate) for summary measures (sW,sA): 
# P_{\bar{g}^*}(sA | sW)/P_{\bar{g}0}(sA | sW)
#-----------------------------------------------------------------------------

# @title Predict h weights under g_0 and g_star using existing m.h.fit model fit
# @name pred.hbars
# @export
# fit models for m_gAi
predict.hbars <- function(newdatnet = NULL, m.h.fit) {
    lbound <- m.h.fit$lbound
    # netA_names <- m.h.fit$m.gAi_vec_g$sA_nms
    # determ_cols_Friend <- m.h.fit$determ_cols_Friend # Should this be saved in m.gAi_vec_g and m.gAi_vec_gstar objects instead?
    # if (!is.null(newdatnet)) {
    #   NetInd_k <- newdatnet$netind_cl$NetInd_k
    #   determ.g_user <- newdatnet$determ.g
    #   determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
    #   determ_cols <- (determ_cols_user | determ_cols_Friend)
    # }
    # use original fitted data for prediction
    if (is.null(newdatnet)) {
      stop("newdatnet argument must be not null; this feature is not implemented")
      # newdatnet <- m.h.fit$cY_mtx_fitted
      # determ_cols <- m.h.fit$determ_cols_fitted
    }
    # PASS ENTIRE newdatnet which will get subset, rather than constructing cY_mtx...
    # if (h_user==FALSE) {
    h_gN <- m.h.fit$summeas.g0$predictAeqa(newdata = newdatnet)
    h_gstar <- m.h.fit$summeas.gstar$predictAeqa(newdata = newdatnet)
    # }
    h_gstar_h_gN <- h_gstar / h_gN
    h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
    h_gstar_h_gN <- bound(h_gstar_h_gN, c(0,1/lbound))
    return(h_gstar_h_gN)
}

# @title Defining and fitting the clever covariate h under g_0 and g_star, i.e. models P(sA[j] | sW,sA[j])
# @name fit.hbars
# @importFrom assertthat assert_that is.count
# @export
# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars <- function(DatNet.ObsP0, est_params_list) {
  # .f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))  # defining the vector of c^A`s that needs evaluation under h(c)
  #---------------------------------------------------------------------------------
  # PARAMETERS FOR LOGISTIC ESTIMATION OF h
  #---------------------------------------------------------------------------------
  O.datnetW <- DatNet.ObsP0$datnetW
  O.datnetA <- DatNet.ObsP0$datnetA
  lbound <- est_params_list$lbound
  max_npwt <- est_params_list$max_npwt # NOT IMPLEMENTED
  ng.MCsims <- est_params_list$ng.MCsims  # replace with p adaptive to k: p <- 100*(2^k)

  sW <- est_params_list$sW
  sA <- est_params_list$sA
  h.g0.sVars <- est_params_list$h.g0.sVars
  h.gstar.sVars <- est_params_list$h.gstar.sVars

  f.gstar <- est_params_list$f.gstar
  f.g0 <- est_params_list$f.g0

  h_g0_SummariesModel <- est_params_list$h_g0_SummariesModel
  if (!is.null(h_g0_SummariesModel)) {
    message("NOTE: Predictions for P(sA|sW) under g0 will be based on the model fit in h_g0_SummariesModel," %+%
            "all modeling settings will be ignored")
  }
  h_gstar_SummariesModel <- est_params_list$h_gstar_SummariesModel
  if (!is.null(h_gstar_SummariesModel)) {
    message("NOTE: Predictions for P(sA^*|sW^*) under f_gstar will be based on the model fit in h_g0_SummariesModel," %+%
      " all modeling settings will be ignored")
  }

  h_logit_sep_k <- est_params_list$h_logit_sep_k # NOT IMPLEMENTED
  # h_user=est_params_list$h_user; h_user_fcn=est_params_list$h_user_fcn; NOT IMPLEMENTED

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sW
  #---------------------------------------------------------------------------------
  # Summary measure names / expressions:
  sW.g0_nms <- h.g0.sVars$predvars
  sW.gstar_nms <- h.gstar.sVars$predvars

  # *****
  # Check that these summary measures exist in O.datnetW$names.sVar
  check.sW.g0.exist <- unlist(lapply(unlist(sW.g0_nms), function(sWname) sWname %in% c(O.datnetW$names.sVar,O.datnetA$names.sVar)))
  if (!all(check.sW.g0.exist)) stop("the following predictors from hform.g0 regression could not be located in sW/sA summary measures: " %+%
                                    paste0(sW.g0_nms[!check.sW.g0.exist], collapse = ","))

  check.sW.gstar.exist <- unlist(lapply(sW.gstar_nms, function(sWname) sWname %in% c(O.datnetW$names.sVar,O.datnetA$names.sVar)))
  if (!all(check.sW.gstar.exist)) stop("the following predictors from hform.gstar regression could not be located in sW/sA summary measures: " %+%
                                    paste0(sW.gstar_nms[!check.sW.gstar.exist], collapse = ","))

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sA
  #---------------------------------------------------------------------------------
  # Summary measure names / expressions:
  sA_nms_g0 <- h.g0.sVars$outvars
  sA_nms_gstar <- h.gstar.sVars$outvars

  # ***********
  # Check that the outcome summary measures defined by h.g0.sVars$outvars and h.gstar.sVars$outvars are equivalent:
  # NOTE: might comment out in the future and allow different summary measures for sA_nms_g0 and sA_nms_gstar.
  # ***********
  for (idx in seq_along(sA_nms_g0)) {
    if (!all(sA_nms_g0[[idx]] %in% sA_nms_gstar[[idx]])) {
      stop("the outcome variable names defined by regressions hform.g0 & hform.gstar are not identical;" %+%
                                            " current implementation requires these to be the same.")
    }
  }

  # ***********
  # Check that these summary measures exist in O.datnetA$names.sVar
  check.sAg0.exist <- unlist(lapply(sA_nms_g0, function(sAname) sAname %in% O.datnetA$names.sVar))
  if (!all(check.sAg0.exist)) stop("the following outcomes from hform.g0 regression could not be located in sA summary measures: " %+%
                                    paste0(sA_nms_g0[!check.sAg0.exist], collapse = ","))

  check.sAgstar.exist <- unlist(lapply(sA_nms_gstar, function(sAname) sAname %in% O.datnetA$names.sVar))
  if (!all(check.sAgstar.exist)) stop("the following outcomes from hform.gstar regression could not be located in sA summary measures: " %+%
                                    paste0(sA_nms_gstar[!check.sAgstar.exist], collapse = ","))


  ##########################################################
  # *********** SUBSETTING ALGORITHM ***********
  ##########################################################
  # NOTE: Subsetting works by var name only (which automatically evaluates as !gvars$misval(var)) for speed & memory efficiency
  # Determine subsetting by !gvars$misval (which observations to include in each univariate regression)
  # (1) For each UNIVARIATE regression (e.g., "A ~ W") specifies a VECTOR of variable names which should be
  #     jointly non-missing in the data, this then defined the observations that go into the design matrix.
  # (2) For each MULTIVARIATE outcome regression with shared predictors (e.g., "A + sumA ~ W"),
  #     this should be a list of length equal to the number of outcomes.
  # (3) For separate regressions with different predictors, i.e., different time-points (e.g., "A_1 + sumA_1 ~ W", "A_2 + sumA_2 ~ W + L_1"),
  #     this should be a list of lists of length equal to the total number of such regressions.
  subsets_expr <- lapply(sA_nms_g0, function(var) lapply(var, function(var) {var}))
  # subsets_expr <- lapply(sA_nms_g0, function(var) {var})

  ##########################################################
  # **** DEPRECATED **** DEFINING SUBSETING EXPRESSIONS (FOR DETERMINISTIC / DEGENERATE sA)
  ##########################################################
  # (1 subset expr per regression P(sA[j]|sA[j-1:0], sW))
  # Old examples of subsetting expressions:
  # based on the variable of gvars$misval (requires passing gvars envir for eval)
  # subset_exprs <- lapply(netvar("determ.g_Friend", c(0:Kmax)), function(var) {var%+%" != "%+%"misval"})
  # based on existing logical determ_g columns (TRUE = degenerate/determ):
  # subset_exprs <- lapply(netvar("determ.g_true", c(0:Kmax)), function(var) {var%+%" != "%+%TRUE})
  #-----------------------------------------------------------

  ##########################################################
  # Summary class params:
  ##########################################
  sA_class <- lapply(sA_nms_g0, function(sA_nms) O.datnetA$type.sVar[sA_nms])
  # sA_class <- O.datnetA$type.sVar[sA_nms_g0[[1]]]

  if (gvars$verbose) {
    message("================================================================")
    message("fitting h_g0 with summary measures: ", "P(" %+% paste(sA_nms_g0, collapse = ",") %+% " | " %+% paste(sW.g0_nms, collapse = ",") %+% ")")
    message("================================================================")
  }

  p_h0 <- ifelse(is.null(f.g0), 1, ng.MCsims)
  if (!is.null(f.g0)) {
    if (gvars$verbose) message("generating DatNet.g0 under known g0")
    DatNet.g0 <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
    DatNet.g0$make.dat.sWsA(p = p_h0, f.g_fun = f.g0, sA.object = sA, DatNet.ObsP0 = DatNet.ObsP0)
    if (gvars$verbose) {
      print("new DatNet.g0$dat.sWsA from known g0: "); print(head(DatNet.g0$dat.sWsA))
    }
  } else {
    DatNet.g0 <- DatNet.ObsP0
  }

  regclass.g0 <- RegressionClass$new(sep_predvars_sets = TRUE,
                                    outvar.class = sA_class,
                                    outvar = sA_nms_g0,
                                    predvars = sW.g0_nms,
                                    subset = subsets_expr)
  regclass.g0$S3class <- "generic"
  # using S3 method dispatch on regclass.g0:
  # time_summeas.g0 <- system.time(
    summeas.g0 <- newsummarymodel(reg = regclass.g0, DatNet.sWsA.g0 = DatNet.g0)
    # )
  # print("time_summeas.g0: "); print(time_summeas.g0)
  
  # summeas.g0 <- SummariesModel$new(reg = regclass.g0, DatNet.sWsA.g0 = DatNet.g0)
  if (!is.null(h_g0_SummariesModel)) {
    # 1) verify h_g0_SummariesModel is consistent with summeas.g0
    assert_that(inherits(h_g0_SummariesModel, "SummariesModel"))
    # 2) deep copy model fits in h_g0_SummariesModel to summeas.g0
    summeas.g0 <- h_g0_SummariesModel$clone(deep=TRUE)
  } else {
    time_summeas.g0.fit <- system.time(
      summeas.g0$fit(data = DatNet.g0)
    )
    print("time_summeas.g0.fit: "); print(time_summeas.g0.fit)
  }


  # *********
  # NEED TO PASS obsdat.sW.sA (observed data sWsA) to predict() funs.
  # If !is.null(f.g_fun) then DatNet.g0$dat.sWsA IS NOT THE OBSERVED data (sWsA), but rather sWsA data sampled under known g_0.
  # Option 1: Wipe out DatNet.g0$dat.sWsA with actually observed data - means that we can't use DatNet.g0$dat.sWsA in the future.
  # Option 2: Create a new class DatNet.Obs of DatNet.sWsA (will be painful)
  # Going with OPTION 1 for now:
  # Already generated DatNet.ObsP0 in tmlenet:
  time_h_gN <- system.time(
    h_gN <- summeas.g0$predictAeqa(newdata = DatNet.ObsP0)
  )
  print("time_h_gN: "); print(time_h_gN)

  if (length(h_gN)!=DatNet.ObsP0$nobs) stop("the IPW weight prediction under g0 return invalid vector length: " %+% length(h_gN))
  # ------------------------------------------------------------------------------------
  # to obtain the decomposed of the above probability by each Anode (marginals):
  # (directly analogous to the g0 component of gmat in ltmle package)
  # gmat.g0 <- matrix(nrow = length(h_gN), ncol = length(summeas.g0$getPsAsW.models()))
  # for (i in seq_along(summeas.g0$getPsAsW.models())) {
  #   gmat.g0[,i] <- summeas.g0$getPsAsW.models()[[i]]$getcumprodAeqa()
  # }
  # cum.gmat.g0 <- matrixStats::rowCumprods(gmat.g0)
  # # message("sum(cum.gmat.g0-h_gN): " %+% sum(cum.gmat.g0[,ncol(cum.gmat.g0)]-h_gN))
  # ------------------------------------------------------------------------------------

  if (gvars$verbose) {
    message("================================================================")
    message("fitting h_gstar based on summary measures: ", "P(" %+% paste(sA_nms_gstar, collapse = ",") %+% " | " %+% paste(sW.gstar_nms, collapse = ",") %+% ")")
    message("================================================================")
  }


  regclass.gstar <- RegressionClass$new(sep_predvars_sets = TRUE,
                                        outvar.class = sA_class,
                                        outvar = sA_nms_gstar,
                                        predvars = sW.gstar_nms,
                                        subset = subsets_expr)
  regclass.gstar$S3class <- "generic"
  # Define Intervals Under g_star to Be The Same as under g0:
  time_summeas.gstar <- system.time(
    summeas.gstar <- newsummarymodel(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.g0)
  )
  print("time_summeas.gstar: "); print(time_summeas.gstar)
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.g0)
  # Define Intervals Under g_star Based on Summary Measures Generated under g_star:
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.gstar)
  # Define Intervals Under g_star Based on Union of Summary Measures under g_star and g0:
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.g0, datnet.gstar = DatNet.gstar)

  DatNet.gstar <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
  time_make.dat.sWsA_gstar <- system.time(
    DatNet.gstar$make.dat.sWsA(p = ng.MCsims, f.g_fun = f.gstar, sA.object = sA, DatNet.ObsP0 = DatNet.ObsP0)
  )
  print("time_make.dat.sWsA_gstar: "); print(time_make.dat.sWsA_gstar)
  if (gvars$verbose) {
    print("Generated new summary measures by sampling A from f_gstar (DatNet.gstar): ")
    print(head(DatNet.gstar$dat.sWsA))
  }

  if (!is.null(h_gstar_SummariesModel)) {
    # 1) verify h_gstar_SummariesModel is consistent with summeas.gstar
    assert_that(inherits(h_gstar_SummariesModel, "SummariesModel"))
    # 2) deep copy the object with model fits to summeas.gstar
    summeas.gstar <- h_gstar_SummariesModel$clone(deep=TRUE)
  } else {
    time_summeas.gstar.fit <- system.time(
      summeas.gstar$fit(data = DatNet.gstar)
    )
    print("time_summeas.gstar.fit: "); print(time_summeas.gstar.fit)
  }

  # All data is now stored in a single global data.table, which is ALWAYS modified IN PLACE
  # For that reason, DatNet.gstar$make.dat.sWsA() call above over-wrote Anodes and sA in Odata$Odata_DT
  # These now need to be restored to the original Anodes and sA under g0:

  # 1) back-up Anodes and sA generated under f.g.star so that we don't have to re-generate them again & then restore old Anodes and sA (under g.0)
  DatNet.gstar$datnetA$Odata$swapAnodes()

 # 2) verify sA's were also restored and if not, regenerate them
  if (!DatNet.ObsP0$datnetA$Odata$restored.sA.Vars)
    DatNet.ObsP0$datnetA$make.sVar(sVar.object = sA)
 
  # browser()
  # DatNet.gstar$datnetA$Odata$OdataDT
  # DatNet.gstar$datnetA$Odata$A_g0_DT
  # DatNet.gstar$datnetA$Odata$sA_g0_DT
 
  time_h_gstar_predict <- system.time(
    h_gstar <- summeas.gstar$predictAeqa(newdata = DatNet.ObsP0)
  )
  print("time_h_gstar_predict: "); print(time_h_gstar_predict)

  if (length(h_gstar)!=DatNet.ObsP0$nobs) stop("the IPW weight prediction under gstar return invalid vector length: " %+% length(h_gstar))

  # ------------------------------------------------------------------------------------
  # to obtain the decomposed probability by each Anode (marginals):
  # (directly analogous to the gstar component of gmat in ltmle package)
  # gmat.gstar <- matrix(nrow = length(h_gstar), ncol = length(summeas.gstar$getPsAsW.models()))
  # for (i in seq_along(summeas.gstar$getPsAsW.models())) {
  #   gmat.gstar[,i] <- summeas.gstar$getPsAsW.models()[[i]]$getcumprodAeqa()
  # }
  # cum.gmat.gstar <- matrixStats::rowCumprods(gmat.gstar)
  # message("sum(gmat.gstar-h_gstar): " %+% sum(cum.gmat.gstar[,ncol(cum.gmat.gstar)]-h_gstar))
  # ------------------------------------------------------------------------------------

  ###########################################
  # 3A) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
  ##########################################
  h_gstar_h_gN <- h_gstar / h_gN
  h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
  h_gstar_h_gN <- bound(h_gstar_h_gN, c(0, 1/lbound))


  ###########################################
  # 3B) Directly evaluating IPTW for static interventions (doesn't need DatNet.gstar)
  ##########################################
  # Ag0mat <- DatNet.ObsP0$get.dat.sWsA(covars = DatNet.ObsP0$nodes$Anodes)
  # Agstarmat <- f.gen.A.star(data = DatNet.ObsP0$dat.sVar, f.g_fun = f.gstar, Anodes = DatNet.ObsP0$nodes$Anodes)
  # h_gstar_I <- Agstarmat == Ag0mat
  # h_gstar_cum <- h_gstar_I[,1]
  # for (i in seq(ncol(h_gstar_I))[-1]) {
  #   h_gstar_cum <- h_gstar_cum*h_gstar_I[,i]
  # }
  # h_gstar_h_gN <- h_gstar_cum / h_gN
  # h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
  # h_gstar_h_gN <- bound(h_gstar_h_gN, c(0, 1/lbound))
  # m.h.fit <- list(summeas.g0 = summeas.g0, summeas.gstar = NULL, lbound = lbound)
  # return(list(h_gstar_h_gN = h_gstar_h_gN, m.h.fit = m.h.fit, DatNet.gstar = NULL))


  ###########################################
  # IPTW BY TIME POINT:
  # NEED TO IMPLEMENT:
  # (1) Apply the bound for each column of the final iptw matrix
  # (2) Be able to easily replace summeas.gstar object with KNOWN probabilities (when gstar is a known user-supplied stochastic intervention)
  # gmat <- cum.gmat.gstar / cum.gmat.g0
  ###########################################
  m.h.fit <- list(summeas.g0 = summeas.g0, summeas.gstar = summeas.gstar, lbound = lbound)
  return(list(h_gstar_h_gN = h_gstar_h_gN, m.h.fit = m.h.fit, DatNet.gstar = DatNet.gstar))
}
