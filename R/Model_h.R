#-----------------------------------------------------------------------------
# Fit and Predict the Efficient IPTW (clever covariate): P_g^*(sA | sW) / P_g0(sA | sW)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# NEW (06/01/15) Get the prob of A^* (known stoch. intervention) from supplied function, fcn_name
#-----------------------------------------------------------------------------
# Get the prob of A^* (known stoch. intervention) from supplied function, fcn_name
# #todo 56 (f.gen.probA.star) +0: 1) rename k to Kmax or remove? 2) remove f_args
f.gen.probA.star <- function(k, df_AllW, fcn_name, f_args = NULL) {
  .f_g_wrapper <- function(k, df_AllW, fcn_name, ...) {
      args0 <- list(k = k, data = df_AllW)
      args <- c(args0, ...)
    do.call(fcn_name, args)
  }
  probA <- .f_g_wrapper(k, df_AllW, fcn_name, f_args)
  return(probA)
}
f.gen.A.star.cont <- function(k, df_AllW, fcn_name, f_args = NULL) {
  f.gen.probA.star(k, df_AllW, fcn_name, f_args)
}

# @title Predict h weights under g_0 and g_star using existing m.h.fit model fit
# @name pred.hbars
# @export
# fit models for m_gAi
pred.hbars <- function(newdatnet = NULL, m.h.fit) {
    NetInd_k <- newdatnet$netind_cl$NetInd_k
    lbound <- m.h.fit$lbound
    # netA_names <- m.h.fit$m.gAi_vec_g$sA_nms
    # determ_cols_Friend <- m.h.fit$determ_cols_Friend # Should this be saved in m.gAi_vec_g and m.gAi_vec_gstar objects instead?
    # if (!is.null(newdatnet)) {
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
    h_gN <- m.h.fit$summeas.g0$predictAeqa(obs.DatNet.sWsA = newdatnet)
    # h_gN <- m.h.fit$summeas.g0$predict(newdata = newdatnet)$predictAeqa(obs.DatNet.sWsA = newdatnet)
    h_gstar <- m.h.fit$summeas.gstar$predictAeqa(obs.DatNet.sWsA = newdatnet)
    # h_gstar <- m.h.fit$summeas.gstar$predict(newdata = newdatnet)$predictAeqa(obs.DatNet.sWsA = newdatnet)
    # }
    h_gstar_gN <- h_gstar / h_gN
    h_gstar_gN[is.nan(h_gstar_gN)] <- 0     # 0/0 detection
    h_gstar_gN <- bound(h_gstar_gN, c(0,1/lbound))
    dat_hest <- data.frame(cY.ID = 0,
                          h_gstar = h_gstar,
                          h_gN = h_gN,
                          h_gstar_gN = h_gstar_gN
                          )
    return(list(dat_hest = dat_hest))
}


# @title Defining and fitting the clever covariate h under g_0 and g_star, i.e. models P(sA[j] | sW,sA[j])
# @name fit.hbars
# @importFrom assertthat assert_that is.count
# @export
# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars <- function(datNetObs, est_params_list) {
  if (gvars$verbose) message("... running fit.hbars() ...")
  .f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))  # defining the vector of c^A's that needs evaluation under h(c)
  #---------------------------------------------------------------------------------
  # PARAMETERS FOR LOGISTIC ESTIMATION OF h
  #---------------------------------------------------------------------------------
  # n <- nrow(data)
  O.datnetW <- datNetObs$datnetW
  O.datnetA <- datNetObs$datnetA

  Kmax <- datNetObs$Kmax
  nodes <- datNetObs$O.datnetW$nodes
  nFnode <- nodes$nFnode
  Anode <- nodes$Anode
  NetInd_k <- datNetObs$O.datnetW$NetInd_k

  lbound <- est_params_list$lbound
  max_npwt <- est_params_list$max_npwt
  ng.MCsims <- est_params_list$ng.MCsims  # replace with p adaptive to k: p <- 100*(2^k)

  sW <- est_params_list$sW
  sA <- est_params_list$sA
  h.sVars <- est_params_list$h.sVars
  h.gstar.sVars <- est_params_list$h.gstar.sVars

  f.gstar <- est_params_list$f.gstar
  f.g0 <- est_params_list$f.g0

  f.g_args <- est_params_list$f.g_args # (HAVE BEEN DISABLED)
  f.g0_args <- est_params_list$f.g0_args # (HAVE BEEN DISABLED)

  h_g0_SummariesModel <- est_params_list$h_g0_SummariesModel
  if (!is.null(h_g0_SummariesModel)) message("h_g0 will not be fit, predicting h_g0 based on the existing fit supplied in h_g0_SummariesModel object")
  h_gstar_SummariesModel <- est_params_list$h_gstar_SummariesModel
  if (!is.null(h_gstar_SummariesModel)) message("h_gstar will not be fit, predicting h_gstar based on the existing fit supplied in h_gstar_SummariesModel object")

  h_logit_sep_k <- est_params_list$h_logit_sep_k # NOT IMPLEMENTED
  # h_user=est_params_list$h_user; h_user_fcn=est_params_list$h_user_fcn; NOT IMPLEMENTED

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sW
  #---------------------------------------------------------------------------------
  # Summary measure names. Will be made specific to g.star or g0.
  sW.g0_nms <- h.sVars$predvars
  sW.gstar_nms <- h.gstar.sVars$predvars
  # *****
  # CHECK THAT THESE SUMMARY MEASURES EXIST IN O.datnetW$names.sVar
  check.sW.g0.exist <- unlist(lapply(sW.g0_nms, function(sWname) sWname %in% O.datnetW$names.sVar))
  assert_that(all(check.sW.g0.exist))

  check.sW.gstar.exist <- unlist(lapply(sW.gstar_nms, function(sWname) sWname %in% O.datnetW$names.sVar))
  assert_that(all(check.sW.gstar.exist))

  if (gvars$verbose) {
    print("check.sW.g0.exist"); print(check.sW.g0.exist)
    print("check.sW.gstar.exist"); print(check.sW.gstar.exist)
  }

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sA
  #---------------------------------------------------------------------------------
  # Summary measure names / expression
  sA_nms <- h.sVars$outvars
  assert_that(all(sA_nms == h.gstar.sVars$outvars))

  # ***********
  # CHECK THAT THESE SUMMARY MEASURES EXIST IN O.datnetW$names.sVar
  check.sA.exist <- unlist(lapply(sA_nms, function(sWname) sWname %in% O.datnetA$names.sVar))
  if (gvars$verbose) {
    print("check.sA.exist"); print(check.sA.exist)
  }
  assert_that(all(check.sA.exist))

  #-----------------------------------------------------------
  # BELOW IS A BUG, all A are assigned to the same bin when trying automatic $detect.sVar.intrvls:
  #-----------------------------------------------------------
  # intvrls <- O.datnetA$detect.sVar.intrvls("A")
  # O.datnetA$set.sVar.intrvls(name.sVar = "A", new.intrvls = intvrls)
  # print("defined intvrls for A: "); print(intvrls)
  # print("TABLE FOR A as an ordinal (all A vals got assigned to cat 2): "); print(table(O.datnetA$discretize.sVar("A")))
  # print("Bins for A: "); print(head(O.datnetA$binirize.sVar("A")))
  # Trying manual intervals for A:
  # O.datnetA$set.sVar.intrvls(name.sVar = "A", new.intrvls = seq(0,1,by=0.1))
  # print(O.datnetA$get.sVar.intrvls("A"))
  # print("Bins for A with manual 10 continuous intervals:")
  # print(table(O.datnetA$discretize.sVar("A")))
  # print(head(O.datnetA$binirize.sVar("A")))
  # print("Stats for sum_1mAW2_nets: ")
  #-----------------------------------------------------------
  # DEFINING SUBSETING EXPRESSIONS (FOR DETERMINISTIC / DEGENERATE sA)
  #-----------------------------------------------------------
  # (1 subset expr per regression P(sA[j]|sA[j-1:0], sW))
  # TO DO: Put this in a separate function (with var as arg + additional args)
  # Old examples of subsetting expressions:
  # subsets_chr <- lapply(sA_nms, function(var) {"!misfun(nFnode)"})  # subsetting by !gvars$misval on sA:
  # subset_exprs <- lapply(sA_nms, function(var) {var%+%" != "%+%"misval"}) # compares to misval constant from gvars envir.    
  # subset_exprs <- lapply(netvar("determ.g_Friend", c(0:Kmax)), function(var) {var%+%" != "%+%"misval"}) # based on the variable of gvars$misval (requires passing gvars envir for eval)
  # subset_exprs <- lapply(netvar("determ.g_true", c(0:Kmax)), function(var) {var%+%" != "%+%TRUE}) # based on existing logical determ_g columns (TRUE = degenerate/determ):
  #-----------------------------------------------------------
  # (1) Capture expression as characeter string: subsetexpr <- deparse(substitute(subsetexpr))
  # (2) Parse the characteer expression into call (make subset expressions into a list of calls (one call per sA[j] in sA))
  # (3) Then substitute the actual var names in the data for generic node names (nFnode, Wnodes, Anode):
  # subsets_chr <- lapply(sA_nms, function(var) {"!misfun("%+%var%+%")"})  # subsetting by !gvars$misval on sA:
  # substitute_list <- lapply(nodes, as.name)
  # subsets_expr <- lapply(subsets_chr, function(subset_chr) {
  #                                         subset_expr <- try(parse(text=subset_chr)[[1]]) # parses chr into a call
  #                                         if(inherits(subset_expr, "try-error")) stop("can't parse the subset formula", call.=FALSE)
  #                                         eval(substitute(substitute(e, env = substitute_list), list(e = subset_expr)))
  #                                       })
  subsets_expr <- lapply(sA_nms, function(var) {var})  # subsetting by !gvars$misval on sA:
  # print("new subsets_expr: "); print(subsets_expr)

  ##########################################
  # Summary class params:
  ##########################################
  sA_class <- O.datnetA$type.sVar[sA_nms]
  # sVartypes <- gvars$sVartypes # <- list(bin = "binary", cat = "categor", cont = "contin")
  # sA_class <- as.list(c(sVartypes$cont, sVartypes$bin, rep_len(sVartypes$bin, length(sA_nms)-2)))

  # sA_class.alt <- lapply(sA_nms, O.datnetA$get.sVar.type)
  # names(sA_class.alt) <- sA_nms
  # print("detected sA_class for O.datnetA$sVar: "); print(sA_class)

  message("================================================================")
  message("fitting h_g0 with summary measures: ", "(" %+% paste(sA_nms, collapse = ",") %+% " | " %+% paste(sW.g0_nms, collapse = ",") %+% ")")
  message("================================================================")

  p_h0 <- ifelse(is.null(f.g0), 1, ng.MCsims)

  if (!is.null(f.g0)) {
    if (gvars$verbose) message("generating DatNet.g0 under known g0")
    DatNet.g0 <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
    DatNet.g0$make.dat.sWsA(p = p_h0, f.g_name = f.g0, f.g_args = f.g0_args, sA.object = sA)
  } else {
    DatNet.g0 <- datNetObs
  }

  regclass.g0 <- RegressionClass$new(outvar.class = sA_class, outvar = sA_nms, predvars = sW.g0_nms, subset = subsets_expr)

  summeas.g0 <- SummariesModel$new(reg = regclass.g0, O.datnetA = DatNet.g0)
  if (!is.null(h_g0_SummariesModel)) {
    message("user supplied model fit for h_g0 is not implemented yet")
    # ...
    # 1) verify h_g0_SummariesModel is consistent with summeas.g0
    # 2) copy model fits in h_g0_SummariesModel to summeas.g0
    # ...
  } else {
    summeas.g0$fit(data = DatNet.g0)
  }

  # *********
  # NEED TO PASS obsdat.sW.sA (observed data sWsA) to predict() funs.
  # If !is.null(f.g_name) then DatNet.g0$dat.sWsA IS NOT THE OBSERVED data (sWsA), but rather sWsA data sampled under known g_0.
  # Option 1: Wipe out DatNet.g0$dat.sWsA with actually observed data - means that we can't use DatNet.g0$dat.sWsA in the future.
  # Option 2: Create a new class DatNet.Obs of DatNet.sWsA - pain in the ass...
  # Going with OPTION 1 for now:
  # Already generated datNetObs in tmlenet:
  h_gN <- summeas.g0$predictAeqa(obs.DatNet.sWsA = datNetObs)
  # *********

  message("================================================================")
  message("fitting h_gstar for summary measures: ", "(" %+% paste(sA_nms, collapse = ",") %+% " | " %+% paste(sW.gstar_nms, collapse = ",") %+% ")")
  message("================================================================")
  DatNet.gstar <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
  DatNet.gstar$make.dat.sWsA(p = ng.MCsims, f.g_name = f.gstar, f.g_args = f.g_args, sA.object = sA)

  if (gvars$verbose) {
    print("DatNet.gstar stored sWsA df: "); print(class(DatNet.gstar$dat.sWsA))
    print(dim(DatNet.gstar$dat.sWsA)); print(head(DatNet.gstar$dat.sWsA));
  }

  # browser()
  # # intervals defined for summary measure sA under g0
  # summeas.g0$getPsAsW.models()[['P(sA|sW).1']]$intrvls
  # # regression class for summary measure sA under g0
  # summeas.g0$getPsAsW.models()[['P(sA|sW).1']]$reg
  # summeas.g0$getPsAsW.models()[['P(sA|sW).1']]$reg$nbins
  # # intervals defined for summary measure net.mean.sA under g0
  # summeas.g0$getPsAsW.models()[['P(sA|sW).2']]$intrvls
  # nbins_gstar <- length(gstar.intrv) - 1
  # bin.nms.sVar_gstar <- DatNet.gstar$bin.nms.sVar("sA", nbins_gstar)
  # print("active bin sVar before calling binirize.sVar: " %+% DatNet.gstar$active.bin.sVar)
  # DatNet.gstar$binirize.sVar(name.sVar = "sA", intervals = gstar.intrv, nbins = nbins_gstar, bin.nms = bin.nms.sVar_gstar)
  # print("active bin sVar after calling binirize.sVar: " %+% DatNet.gstar$active.bin.sVar)
  # print(head(DatNet.gstar$dat.sVar, 5))
  # print(head(cbind(DatNet.gstar$ord.sVar, DatNet.gstar$dat.bin.sVar), 5))
  # print("freq count for transformed ord.sVar: "); print(table(DatNet.gstar$ord.sVar))

  regclass.gstar <- RegressionClass$new(outvar.class = sA_class,
                                        outvar = sA_nms,
                                        predvars = sW.gstar_nms,
                                        subset = subsets_expr
                                        # max_nperbin = as.integer(getopt("maxNperBin")*ng.MCsims)
                                        )
  # Define Intervals Under g_star to Be The Same as under g0:
  summeas.gstar <- SummariesModel$new(reg = regclass.gstar, O.datnetA = DatNet.g0)
  # Define Intervals Under g_star Based on Summary Measures Generated under g_star:
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, O.datnetA = DatNet.gstar)
  # Define Intervals Under g_star Based on Union of Summary Measures under g_star and g0:
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, O.datnetA = DatNet.g0, datnet.gstar = DatNet.gstar)

  if (!is.null(h_gstar_SummariesModel)) {
    message("user supplied model fit for h_gstar is not implemented yet")
    # ...
    # 1) verify h_gstar_SummariesModel is consistent with summeas.gstar
    # 2) copy model fits in h_gstar_SummariesModel to summeas.gstar
    # ...
  } else {
    summeas.gstar$fit(data = DatNet.gstar)
  }

  # summeas.gstar$predict(newdata = datNetObs)
  h_gstar <- summeas.gstar$predictAeqa(obs.DatNet.sWsA = datNetObs)
  # # new_h_gstar <- h_gstar
  # # newint_pred <- as.vector(new_h_gstar)
  # # old_h_gstar <- h_gstar
  # # oldint_pred <- as.vector(old_h_gstar)
  
  # # print(head(datNetObs$dat.sWsA))
  # dat_pred <- cbind(datNetObs$dat.sWsA[,c("sA")], oldint_pred=oldint_pred, newint_pred=newint_pred, diff=newint_pred-oldint_pred)
  # # head(dat_pred, 1000)
  # # dat_pred[42,]
  # # dat_pred[11172,]
  # (g0_intrvls <- summeas.g0$getPsAsW.models()[['P(sA|sW).1']]$intrvls)
  # (gstar_intrvls <- summeas.gstar$getPsAsW.models()[['P(sA|sW).1']]$intrvls)
  # # (combine_int <- sort(union(g0_intrvls, gstar_intrvls)))

  # 1/diff(g0_intrvls)
  # 1/diff(gstar_intrvls)
  # # datNetObs$mat.bin.sVar[11172,]
  # # datNetObs$ord.sVar[11172]
  # # h_gN[11172]
  # # # obs.DatNet.sWsA$get.sVar.bw(name.sVar = self$outvar, intervals = self$intrvls)
  # # print(mean(h_gstar_gN))
  # # which(is.nan(h_gstar_gN))
  # # which(h_gstar_gN > 10000)
  # # dat_pred[which(abs(new_wts-old_wts) > 0.4)[17],]
  # # new_wts[which(abs(new_wts-old_wts) > 0.4)[17]]
  # # old_wts[which(abs(new_wts-old_wts) > 0.4)[17]]
  # # h_gstar[which(abs(new_wts-old_wts) > 0.4)[17]]
  # # new_h_gstar[which(abs(new_wts-old_wts) > 0.4)[17]]
  # # old_h_gstar[which(abs(new_wts-old_wts) > 0.4)[17]]
  # # h_gN[which(abs(new_wts-old_wts) > 0.4)[17]]
  # # mean(new_h_gstar)
  # # mean(old_h_gstar)

  # mean(oldint_pred[selidx]/h_gN[selidx])
  # mean(newint_pred[selidx]/h_gN[selidx])

  # selidx <- datNetObs$dat.sWsA[,"W1"]==1 & datNetObs$dat.sWsA[,"W2"]==1 &  datNetObs$dat.sWsA[,"W3"]==1
  # sum(selidx)
  # length(bigdiffidx)
  # dat_pred[bigdiffidx & selidx,][1:20,]

  # mean(datNetObs$dat.sWsA[selidx,c("sA")])
  # plot(density(DatNet.gstar$dat.sWsA[selidx,c("sA")]), col="red")
  # lines(density(datNetObs$dat.sWsA[selidx,c("sA")]))


  ###########################################
  # 3) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
  ##########################################
  h_gstar_gN <- h_gstar / h_gN
  h_gstar_gN[is.nan(h_gstar_gN)] <- 0     # 0/0 detection
  h_gstar_gN <- bound(h_gstar_gN, c(0, 1/lbound))

  # # new_wts <- h_gstar_gN
  # # mean(new_wts)
  # old_wts <- h_gstar_gN
  # mean(old_wts)
  # mean(old_wts); mean(new_wts)
  # bigdiffidx <- abs(new_wts-old_wts) > 0.4
  # table(names(which(bigdiffidx)))
  # # clearly most of disagreement is occuring at the rightmost interval end:
  # # 100%  50%  85%  90%  95%
  # #  996   62   55  190  439


  dat_hest <- data.frame(cY.ID = .f.mkstrNet(datNetObs$dat.sWsA),
                        h_gstar = h_gstar,
                        h_gN = h_gN,
                        h_gstar_gN = h_gstar_gN
                        )

  m.h.fit <- list(summeas.g0 = summeas.g0,
                summeas.gstar = summeas.gstar,
                lbound=lbound
                )

  return(list(dat_hest = dat_hest, m.h.fit = m.h.fit, datNetObs = datNetObs, DatNet.gstar = DatNet.gstar))
}
