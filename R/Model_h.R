gvars <- new.env(parent = emptyenv())
gvars$verbose <- FALSE      # verbose mode (print all messages)
gvars$misval <- NA_integer_ # the default missing value for observations
# gvars$misval <- 1L
# gvars$misval <- -.Machine$integer.max
gvars$misXreplace <- 0L     # the default replacement value for misval that appear in the design matrix
gvars$tolerr <- 10^-12      # tolerance error: assume for abs(a-b) < gvars$tolerr => a = b
gvars$maxncats <- 5L        # max number of categories for sA[j]. If more, automatically considered continuous
# gvars$maxncats <- 10L
gvars$nbins <- 10L          # default n bins for continous var
gvars$sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")

# returns a function (alternatively a call) that tests for missing values in sA / sW
testmisfun <- function() {
  if (is.na(gvars$misval)) {
    return(is.na)
  } else if (is.null(gvars$misval)){
    return(is.null)
  } else if (is.integer(gvars$misval)) {
    return(function(x) {x==gvars$misval})
  } else {
    return(function(x) {x%in%gvars$misval})
  }
}
get.misval <- function() {
  gvars$misfun <- testmisfun()
  gvars$misval
}
set.misval <- function(gvars, newmisval) {
  oldmisval <- gvars$misval
  gvars$misval <- newmisval
  gvars$misfun <- testmisfun()    # EVERYTIME gvars$misval HAS CHANGED THIS NEEDS TO BE RESET/RERUN.
  invisible(oldmisval)
}
set.nbins <- function(nbins) {
  old.nbins <- gvars$nbins
  gvars$nbins <- nbins
  invisible(old.nbins)
}
set.maxncats <- function(maxncats) {
  old.maxncats <- gvars$maxncats
  gvars$maxncats <- maxncats
  invisible(old.maxncats)
}
gvars$misfun <- testmisfun()

#-----------------------------------------------------------------------------
# Class Membership Tests
#-----------------------------------------------------------------------------
is.DatNet.sWsA <- function(DatNet.sWsA) "DatNet.sWsA"%in%class(DatNet.sWsA)
is.DatNet <- function(DatNet) "DatNet"%in%class(DatNet)

#-----------------------------------------------------------------------------
# ALL NETWORK VARIABLE NAMES MUST BE CONSTRUCTED BY CALLING THIS FUNCTION.
# In the future might return the network variable (column vector) itself.
# Helper function that for given variable name (varnm) and friend index (fidx) 
# returns the characeter name of that network variable varnm[fidx], 
# for fidx = 0 (var itself), ..., kmax. fidx can be a vector, in which case a 
# character vector of network names is returned. If varnm is also a vector, a 
# character vector for all possible combinations of (varnm x fidx) is returned.
#-----------------------------------------------------------------------------
# OUTPUT format: Varnm_net.j:
netvar <- function(varnm, fidx) {
  cstr <- function(varnm, fidx) {
    slen <- length(fidx)
    rstr <- vector(mode = "character", length = slen)
    netidxstr <- ! (fidx %in% 0L)
    rstr[netidxstr] <- stringr::str_c('_netF', fidx[netidxstr])  # vs. 1
    # rstr[netidxstr] <- str_c('.net.', fidx[netidxstr])  # vs. 2
    return(stringr::str_c(varnm, rstr))
  }
  if (length(varnm) > 1) {
    return(unlist(lapply(varnm, cstr, fidx)))
  } else {
    return(cstr(varnm, fidx))
  }
}
# Examples:
# netvar("A", (0:5))
# netvar("A", c(0:5, 0, 3))
# netvar(c("A", "W"), c(0:5, 0, 3))
# netvar(c("A", "W"), c(0:5, 0, 3))


#-----------------------------------------------------------------------------
# Bound g(A|W) probability within supplied bounds
#-----------------------------------------------------------------------------
bound <- function(x, bounds){
	x[x<min(bounds)] <- min(bounds)
	x[x>max(bounds)] <- max(bounds)
	return(x)
}

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

f.gen.A.star <- function(k, df_AllW, fcn_name, f_args = NULL) {
  n <- nrow(df_AllW)
  rbinom(n, 1, f.gen.probA.star(k, df_AllW, fcn_name, f_args))
}

# (DEPRECATED) NEEDS TO BE REPLACED WITH BinOutModel based prediction
#-----------------------------------------------------------------------------
# get the prob of A (under g_0) using fit g_N, estimated regression model for g0;
#-----------------------------------------------------------------------------
.f.gen.probA_N <- function(df, deterministic, m.gN) {
    SuppressGivenWarnings({
    	g_N <- predict(m.gN, newdata=df, type="response")
					},
					GetWarningsToSuppress())
    #************************************************
  	g_N[deterministic] <- 0
    #************************************************
  	return(g_N)
}

#------------------------------------------------------------------------------
# NOT FINISHED
# # todo 65 (iptw_est) +0: rework iptw_est functions for new package structure
# IPTW ESTIMATOR (est Y_g_star based on weights g_star(A|W)/g_N(A|W) )
#------------------------------------------------------------------------------
iptw_est <- function(k, data, node_l, m.gN, f.gstar, f.g_args, family="binomial", NetInd_k, lbound=0, max_npwt=50, f.g0=NULL, f.g0_args=NULL, iidIPTW=FALSE) {
	n <- nrow(data)
	netW <- NULL
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
	for (Wnode in node_l$Wnodes) {
	 	netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[, Wnode], Wnode)))
	}
  cA.mtx <- cbind(netW, subset(data, select=nFnode))
	indA <- data.frame(.f.allCovars(k, NetInd_k, data[, Anode], Anode))

  determ.g <- data$determ.g
  determ.g_vals <- data[determ.g, Anode]
	# predict g*(A=1|W):
	pA_gstar <- f.gen.probA.star(k, cA.mtx, f.gstar, f.g_args)
	netpA_gstar <- .f.allCovars(k, NetInd_k, pA_gstar, Anode)

	# calculate likelihoods P_*(A=a|W):
  if (iidIPTW) { # for iid IPTW, use only A_i, no A_j, for j\in F_i:
    indA <- indA[, Anode, drop=FALSE]
    netpA_gstar <- netpA_gstar[, Anode, drop=FALSE]
  }

  gstar_A <- .f.cumprod.matrix(indA, netpA_gstar) # calculate likelihoods P_0(A=a|W) under g_star:
  #************************************************
  if (is.null(f.g0)) { #If g_0 is unknown, use logistic fit m.gN
    # print("RUNNING IPTW ON FITTED g0_N"); print(m.gN)
    pA_g0N <- .f.gen.probA_N(cA.mtx, determ.g, m.gN)
  }
  else {   # If g_0 is known get true P_g0
    # print("RUNNING IPTW ON TRUE g0")
    pA_g0N <- f.gen.probA.star(k, cA.mtx, f.g0, f.g0_args)
  }
  #************************************************
	netpA_g0 <- .f.allCovars(k, NetInd_k, pA_g0N, Anode)

  # for iid IPTW, use only A_i, no A_j, for j\in F_i:
  if (iidIPTW) {
    indA <- indA[,Anode,drop=FALSE]
    netpA_g0 <- netpA_g0[,Anode,drop=FALSE]
  }
  # print("head(indA)"); print(head(indA))
  # print("head(netpA_g0)"); print(head(netpA_g0))

  # calculate likelihoods P_0(A=a|W):
  g0_A <- .f.cumprod.matrix(indA, netpA_g0)
  ipweights <- gstar_A / g0_A
  ipweights[is.nan(ipweights)] <- 0     # 0/0 detection
  # ipweights <- bound(ipweights, c(0,10^6))
  # lower bound g0 by lbound
  ipweights <- bound(ipweights, c(0,1/lbound))

  # scale weights by total contribution (instead of bounding by consts):
  # cap the prop weights scaled at max_npwt (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000)
  # ipweights <- scale_bound(ipweights, max_npwt, n)
  # print("iptw wts range after bound"); print(summary(ipweights))
	return(ipweights)
}

#' @title Predict the h weights under g_0 and g_star using existing m.h.fit model fit
#' @docType function
#' @name pred.hbars
##' @export
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
    if (is.null(newdatnet)) {    # use original fitted data for prediction
      # ... NOT IMPLEMENTED ...
      # newdatnet <- m.h.fit$cY_mtx_fitted
      # determ_cols <- m.h.fit$determ_cols_fitted
    }
    # PASS ENTIRE newdatnet which will get subset, rather than constructing cY_mtx...
    # if (h_user==FALSE) {
      h_gN <- m.h.fit$summeas.g0$predict(newdata = newdatnet)$predictAeqa(obs.DatNet.sWsA = newdatnet)
      h_gstar <- m.h.fit$summeas.gstar$predict(newdata = newdatnet)$predictAeqa(obs.DatNet.sWsA = newdatnet)
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


#' @title Defining and fitting the clever covariate h under g_0 and g_star, i.e. models P(sA[j]|sW,\bar{sA}[j])
#' @docType function
#' @name fit.hbars
# @importFrom assertthat assert_that is.count
##' @export
# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars <- function(datNetObs, est_params_list) {
  message("... running fit.hbars() ...")  
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
  h.sVars = est_params_list$h.sVars
  h.gstar.sVars = est_params_list$h.gstar.sVars

  f.gstar <- est_params_list$f.gstar
  f.g0 <- est_params_list$f.g0

  f.g_args <- est_params_list$f.g_args # TO BE REmoVED
  f.g0_args <- est_params_list$f.g0_args # TO BE REmoVED
  h_logit_sep_k <- est_params_list$h_logit_sep_k # NOT IMPLEMENTED
  # h_user=est_params_list$h_user; h_user_fcn=est_params_list$h_user_fcn; NOT IMPLEMENTED

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sW
  #---------------------------------------------------------------------------------
  # Summary measure names. Will be made specific to g.star or g0.
  sW.g0_nms <- h.sVars$predvars
  sW.gstar_nms <- h.gstar.sVars$predvars
  # sW_nms <- O.datnetW$names.sVar

  # *****
  # CHECK THAT THESE SUMMARY MEASURES EXIST IN O.datnetW$names.sVar
  check.sW.g0.exist <- unlist(lapply(sW.g0_nms, function(sWname) sWname %in% O.datnetW$names.sVar))
  print("check.sW.g0.exist"); print(check.sW.g0.exist)
  assert_that(all(check.sW.g0.exist))

  check.sW.gstar.exist <- unlist(lapply(sW.gstar_nms, function(sWname) sWname %in% O.datnetW$names.sVar))
  print("check.sW.gstar.exist"); print(check.sW.gstar.exist)
  assert_that(all(check.sW.gstar.exist))
  # *****

    # print("type.sVar: "); str(O.datnetW$type.sVar)
    # print("names.c.sVar: "); print(O.datnetW$names.c.sVar)
  # Detect intervals for continous covars (Don't really need to do that for netW)
    # O.datnetW$def_cbin_intrvls()
    # print("Detected intervals: "); print(O.datnetW$cbin_intrvls)
    # print("Detected nbins: "); print(O.datnetW$all.nbins)
  #---------------------------------------------------------------------------------
  # Getting OBSERVED sA (both, h.sVars or h.gstar.sVars, have to have the same sAs!)
  #---------------------------------------------------------------------------------
  # Summary measure names / expression
  sA_nms <- h.sVars$outvars


  assert_that(all(sA_nms == h.gstar.sVars$outvars))
  # *****
  # CHECK THAT THESE SUMMARY MEASURES EXIST IN O.datnetW$names.sVar
  check.sA.exist <- unlist(lapply(sA_nms, function(sWname) sWname %in% O.datnetA$names.sVar))
  print("check.sA.exist"); print(check.sA.exist)
  assert_that(all(check.sA.exist))
  # *****

  # sA_nms <- O.datnetA$names.sVar
  # ***********
  # Detect intervals for continous covars in sVar
    # O.datnetA$def_cbin_intrvls()
    # print("Detected types: "); str(O.datnetA$type.sVar)
    # print("Detected intervals: "); print(O.datnetA$cbin_intrvls)
    # print("Detected nbins: "); print(O.datnetA$all.nbins)

  #-----------------------------------------------------------
  # Turn A into cont type and add some intervals to it
  #-----------------------------------------------------------
  message("Testing manually to set A as cont:")

  # O.datnetA$set.sVar.type(name.sVar = "A", new.type = "contin")
  message("This is a bug, all A are assigned to the same bin when trying automatic $detect.sVar.intrvls:")
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

  print("Stats for sum_1mAW2_nets: ")
  # print(O.datnetA$get.sVar.type("sum_1mAW2_nets"))
  # print(O.datnetA$get.sVar.intrvls("sum_1mAW2_nets"))
  # print(O.datnetA$nbins.sVar("sum_1mAW2_nets"))
  # print(O.datnetA$bin.nms.sVar("sum_1mAW2_nets"))
  # print(table(O.datnetA$discretize.sVar("sum_1mAW2_nets")))
  # print(head(O.datnetA$binirize.sVar("sum_1mAW2_nets")))
  # print("detect sVar-specific intervals: sum_1mAW2_nets"); O.datnetA$detect.sVar.intrvls("sum_1mAW2_nets")
  # TO SET sVar intervals for specific sVar name:
  # O.datnetA$set.sVar.intrvls("sum_1mAW2_nets", O.datnetA$detect.sVar.intrvls("sum_1mAW2_nets"))
  # print("currently defined intervals: "); print(O.datnetA$cbin_intrvls)
  # print("detect and set ALL intervals: "); O.datnetA$def_cbin_intrvls()
  # print("intervals after automatic detection: "); print(O.datnetA$cbin_intrvls)
  # browser()

  #-----------------------------------------------------------
  # DEFINING SUBSETING EXPRESSIONS (FOR DETERMINISTIC / DEGENERATE sA)
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
  # TO DO: Replace var with the index (column) such as, which(colnames(dat)%in%var)
  # Then can use evaluator even when dat is matrix
  subsets_chr <- lapply(sA_nms, function(var) {"!misfun("%+%var%+%")"})  # subsetting by !gvars$misval on sA:
  substitute_list <- lapply(nodes, as.name)
  subsets_expr <- lapply(subsets_chr, function(subset_chr) {
                                          subset_expr <- try(parse(text=subset_chr)[[1]]) # parses chr into a call
                                          if(inherits(subset_expr, "try-error")) stop("can't parse the subset formula", call.=FALSE)
                                          eval(substitute(substitute(e, env = substitute_list), list(e = subset_expr)))
                                        })
  # print("subsets_expr: "); print(subsets_expr)

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
  message("fitting h_g0 with summary measures: ")
  message("(" %+% paste(sA_nms, collapse = ",") %+% " | " %+% paste(sW.g0_nms, collapse = ",") %+% ")")
  message("================================================================")

  p_h0 <- ifelse(is.null(f.g0), 1, ng.MCsims)
  
# ******************************************************************************************************************
#                         **************** NOTE. VERY IMPORTANT ****************
# Note the order in which SummariesModel$new & DatNet.g0$make.dat.sWsA are called is very important!!!!
# If relyng on automatic interval detection for inside ContinSummaryModel$new (O.datnetA$set.sVar.intrvls & O.datnetA$detect...)
# for cont sVar (that hasn't had its intervals defined yet), calling DatNet.g0$make.dat.sWsA BEFORE SummariesModel$new is bad!
# Will result in DatNet.g0 copying & using old interval definitions from O.datnetA,
# even though new interval defns will be later saved inside O.datnetA, DatNet.g0 will not see them.
# When binirize on DatNet.g0 is called from summeas.g0$fit it will use whatever intervals were copied at time of $make.dat.sWsA
# This is due to: 
# (1) side-effects based programming (interval defn in SummariesModel$new is a side-effect)
# (2) $binirize() being part of DatNet class and not DatNet.sWsA (and hence not seeing intervals in O.datnetA)
# (3) Another related problem is that datNetObs$make.dat.sWsA is called in tmlenet(), where no intervals were defined
  # as a result DatNet.g0 below has NO INTERVALS SAVED => Now forcing $copy.cbin.intrvls() on it which fixes the issue.
# TO FIX THIS:
# 1) Call data$copy.cbin.intrvls() from ContinSummaryModel$fit => MIGHT BE THE EASIEST FIX (THIS IS WHAT IS CURRENTLY DONE)
# 2) Move $binirize functions to DatNet.sWsA and make binirize always first check and copy interval defns from O.datnetA
# 3) Define ALL intervals AT THE SAME TIME as creating O.datnetA$new (possibly redefining when creating datnetA under known g0).
    # Might not fully fix the issue, since DatNet.sWsA is created inside both, tmlenet and fit.h functions
# ******************************************************************************************************************

  datNetObs$copy.cbin.intrvls()
  if (!is.null(f.g0)) {
    message("generating DatNet.g0 under known g0")
    # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
    DatNet.g0 <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
    DatNet.g0$make.dat.sWsA(p = p_h0, f.g_name = f.g0, f.g_args = f.g0_args, sA.object = sA)
  } else {
    DatNet.g0 <- datNetObs
  }

  # (1) TO BE REPLACED WITH RegressionClass object that defines sA_class, sA_nms, sW_nms
  # (2) #todo 40 (SummariesModel$new) +0: PASS datNetObs to constructor INSTEAD OF O.datnetA?
  summeas.g0 <- SummariesModel$new(sA_class = sA_class,
                                    sA_nms = sA_nms,
                                    sW_nms = sW.g0_nms,
                                    subset = subsets_expr,
                                    O.datnetA = O.datnetA)

  # print("all intervals DatNet.g0$datnetA: "); print(DatNet.g0$datnetA$cbin_intrvls)
  # print("all intervals (for O.datnetA & O.datnetW) in DatNet.g0: "); print(DatNet.g0$cbin_intrvls)
  # browser()
  summeas.g0$fit(data = DatNet.g0)

  # *********
  # NEED TO PASS obsdat.sW.sA (observed data sWsA) to predict() funs.
  # predict() expects obsdat.sW.sA to be also of class DatNet.sWsA and will call on the same methods of DatNet.sWsA as does fit() (constructing bins, asking for data.frame, etc)
  # If !is.null(f.g_name) then DatNet.g0$dat.sWsA IS NOT THE OBSERVED data (sWsA), but rather sWsA data sampled under known g_0.
  # Option 1: Wipe out DatNet.g0$dat.sWsA with actually observed data - means that we can't use DatNet.g0$dat.sWsA in the future.
  # Option 2: Create a new class DatNet.Obs of DatNet.sWsA - pain in the ass...
  # Going with OPTION 1 for now:

  # Already generated datNetObs in tmlenet:
  summeas.g0$predict(newdata = datNetObs) # DOESN'T HAVE TO BE CALLED if (is.null(f.g0)), since PREDICATIONS ARE ALREADY SAVED for obsdat.sW.sA:
  # *********
  h_gN <- summeas.g0$predictAeqa(obs.DatNet.sWsA = datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  # *********

  message("================================================================")
  message("fitting h_gstar for summary measures: ")
  message("(" %+% paste(sA_nms, collapse = ",") %+% " | " %+% paste(sW.gstar_nms, collapse = ",") %+% ")")
  message("================================================================")

  # temp setting: ng.MCsims <- 10
  DatNet.gstar <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
  DatNet.gstar$make.dat.sWsA(p = ng.MCsims, f.g_name = f.gstar, f.g_args = f.g_args, sA.object = sA)

  print("DatNet.gstar stored sWsA df: "); print(class(DatNet.gstar$dat.sWsA))
  print(dim(DatNet.gstar$dat.sWsA)); print(head(DatNet.gstar$dat.sWsA));
  # print("all intervals (for O.datnetA & O.datnetW) in DatNet.gstar: "); print(DatNet.gstar$cbin_intrvls)
  # print("all var types (for O.datnetA & O.datnetW) in DatNet.gstar: "); str(DatNet.gstar$type.sVar)
  # DatNet.gstar$names.sVar
  # head(DatNet.gstar$dat.sVar)

  summeas.gstar <- SummariesModel$new(sA_class = sA_class,
                                      sA_nms = sA_nms,
                                      sW_nms = sW.gstar_nms,
                                      subset = subsets_expr,
                                      O.datnetA = O.datnetA)

  summeas.gstar$fit(data = DatNet.gstar)
  summeas.gstar$predict(newdata = datNetObs)
  h_gstar <- summeas.gstar$predictAeqa(obs.DatNet.sWsA = datNetObs)

  ###########################################
  # 3) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
  ##########################################
  h_gstar_gN <- h_gstar / h_gN
  h_gstar_gN[is.nan(h_gstar_gN)] <- 0     # 0/0 detection
  h_gstar_gN <- bound(h_gstar_gN, c(0,1/lbound))

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
