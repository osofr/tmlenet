library(stringr)

`%+%` <- function(a, b) paste0(a, b)  # cat function
gvars <- new.env(parent = emptyenv())

# gvars$misval <- 1L
# gvars$misval <- -.Machine$integer.max
gvars$misval <- NA_integer_
gvars$misXreplace <- 0L
testmisfun <- function() {  # returns a function (alternatively a call) that tests for missing values in sA / sW
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
gvars$misfun <- testmisfun()
setmisval <- function(gvars, newmisval) {
  oldmisval <- gvars$misval
  gvars$misval <- newmisval
  gvars$misfun <- testmisfun()    # EVERYTIME gvars$misval HAS CHANGED THIS NEEDS TO BE RESET/RERUN.
  invisible(oldmisval)
}
getmisval <- function() {
  gvars$misfun <- testmisfun()
  gvars$misval
}


#-----------------------------------------------------------------------------
# ALL NETWORK VARIABLE NAMES MUST BE CONSTRUCTED BY CALLING THIS FUNCTION.
# In the future might return the network variable (column vector) itself.
# Helper function that for given variable name (varnm) and friend index (fidx) 
# returns the characeter name of that network variable varnm[fidx], 
# for fidx = 0 (var itself), ..., kmax. fidx can be a vector, in which case a 
# character vector of network names is returned. If varnm is also a vector, a 
# character vector for all possible combinations of (varnm x fidx) is returned.
#-----------------------------------------------------------------------------
netvar <- function(varnm, fidx) { # OUTPUT format: netVarnm_j
  cstr <- function(varnm, fidx) {
    slen <- length(fidx)
    lstr <- vector(mode = "character", length = slen)
    rstr <- vector(mode = "character", length = slen)
    netidxstr <- !(fidx %in% 0L)
    lstr[netidxstr] <- "net"    
    rstr[netidxstr] <- str_c('_', fidx[netidxstr])
    return(str_c(lstr, varnm, rstr))
  }
  if (length(varnm) > 1) {
    return(unlist(lapply(varnm, cstr, fidx)))
  } else {
    return(cstr(varnm, fidx))
  }
}
netvar2 <- function(varnm, fidx) { # OUTPUT format: Varnm_net.j
  cstr <- function(varnm, fidx) {
    slen <- length(fidx)
    rstr <- vector(mode = "character", length = slen)
    netidxstr <- !(fidx %in% 0L)
    rstr[netidxstr] <- str_c('_netF', fidx[netidxstr])  # vs. 1
    # rstr[netidxstr] <- str_c('.net.', fidx[netidxstr])  # vs. 2
    return(str_c(varnm, rstr))
  }
  if (length(varnm) > 1) {
    return(unlist(lapply(varnm, cstr, fidx)))
  } else {
    return(cstr(varnm, fidx))
  }
}
# Examples:
# netvar("A", (0:5))
# netvar2("A", c(0:5, 0, 3))
# netvar(c("A", "W"), c(0:5, 0, 3))
# netvar2(c("A", "W"), c(0:5, 0, 3))


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
f.gen.probA.star <- function(k, df_AllW, fcn_name, f_args=NULL) { # Get the prob of A^* (known stoch. intervention) from supplied function, fcn_name
  .f_g_wrapper <- function(k, df_AllW, fcn_name, ...) {
      args0 <- list(k=k, data=df_AllW)
      args <- c(args0, ...)
    do.call(fcn_name, args)
  }
  probA <- .f_g_wrapper(k, df_AllW, fcn_name, f_args)
  return(probA)
}
f.gen.A.star <- function(k, df_AllW, fcn_name, f_args=NULL) {
  n <- nrow(df_AllW)
  rbinom(n, 1, f.gen.probA.star(k, df_AllW, fcn_name, f_args))
}

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
.f.gen.A_N <- function(df, deterministic, m.gN) { # sample treatments with probA from above fcn
	n <- nrow(df)
  	rbinom(n, 1, .f.gen.probA_N(df, deterministic, m.gN))
}

#-----------------------------------------------------------------------------
# return entire network matrix from indiv. covar (Var) + covariate itself as first column
#-----------------------------------------------------------------------------
.f.allCovars <- function(k, NetInd_k, Var, VarNm, misval = 0L) {
	n <- length(Var) 
	NetInd_k <- matrix(NetInd_k, nrow=n, ncol=k)
	netVar_names <- NULL
	netVar_full <- NULL
	d <- matrix(0L, nrow=n, ncol = k+1)
	d[ , 1] <- Var
	d[ , c(2:(k+1))] <- apply(NetInd_k, 2, function(k_indx) {
                    											netVar <- Var[k_indx]
                                          netVar[is.na(netVar)] <- misval
                    											return(netVar)
                    											})
  if (k>0) netVar_names <- netvar2(varnm = VarNm, fidx = c(1:k))
	Var_names <- c(VarNm, netVar_names)
	colnames(d) <- Var_names
	return(d)
}

#-----------------------------------------------------------------------------
# Fit glm based on formula in "form"
#-----------------------------------------------------------------------------
.f.est <- function(d, form, family) {
  ctrl <- glm.control(trace=FALSE, maxit=1000)
    SuppressGivenWarnings({
              m <- glm(as.formula(form), 
                  data=d, 
                  family=family, 
                  control=ctrl)
              }, GetWarningsToSuppress())
    return(m)
}
#------------------------------------------------------------------------------
# IPTW ESTIMATOR (est Y_g_star based on weights g_star(A|W)/g_N(A|W) )
#------------------------------------------------------------------------------
iptw_est <- function(k, data, node_l, m.gN, f.g.star, f.g_args, family="binomial", NetInd_k, lbound=0, max_npwt=50, f.g0=NULL, f.g0_args=NULL, iidIPTW=FALSE) {	
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
	pA_gstar <- f.gen.probA.star(k, cA.mtx, f.g.star, f.g_args)
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

pred.hbars.new <- function(new_data=NULL, fit_h_reg_obj, NetInd_k) {
    lbound <- fit_h_reg_obj$lbound
    netA_names <- fit_h_reg_obj$m.gAi_vec_g$sA_nms

    determ_cols_Friend <- fit_h_reg_obj$determ_cols_Friend # Should this be saved in m.gAi_vec_g and m.gAi_vec_gstar objects instead?

    if (!is.null(new_data)) {
      determ.g_user <- new_data$determ.g
      determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
      determ_cols <- (determ_cols_user | determ_cols_Friend)
    }
    if (is.null(new_data)) {    # use original fitted data for prediction
      # ... NOT IMPLEMENTED ... 
      # new_data <- fit_h_reg_obj$cY_mtx_fitted
      # determ_cols <- fit_h_reg_obj$determ_cols_fitted
    }
    indA <- as.matrix(new_data[, netA_names])
    # PASS ENTIRE new_data which will get subset, rather than constructing cY_mtx...
    new_data <- cbind(determ_cols, new_data)

    print("head(new_data):"); print(head(new_data))
    # cY_mtx <- cbind(determ_cols, as.matrix(new_data[, c(node_l$nFnode, W_nms, netA_names)]))

    if (h_user==FALSE) {
      h_vec.g0.new <- fit_h_reg_obj$m.gAi_vec_g$predict(newdata = new_data)$predictAeqa(indA = indA)$getcumprodAeqa()
      h_vec.gstar.new <- fit_h_reg_obj$m.gAi_vec_gstar$predict(newdata = new_data)$predictAeqa(indA = indA)$getcumprodAeqa()
      #---------------------------------------------------------------------
    }
    h_tilde.new <- h_vec.gstar.new / h_vec.g0.new
    h_tilde.new[is.nan(h_tilde.new)] <- 0     # 0/0 detection
    h_tilde.new <- bound(h_tilde.new, c(0,1/lbound))
    df_h_bar_vals <- data.frame(cY.ID = 0, 
                                h.star_c = h_vec.gstar.new, 
                                h_c = h_vec.g0.new,
                                h = h_tilde.new
                                )
    return(list(df_h_bar_vals=df_h_bar_vals))
}

# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars.new <- function(data, h_fit_params) {    
    .f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" ")) # defining the vector of c^A's that needs evaluation under h(c) 
    #---------------------------------------------------------------------------------
    # GENERATE SAMPLES OF (dat.sW, dat.sA), where sA is sampled under f.g_name, of size p*nobs, for p>=1.
    # Returns observed (sA,sW) data if is.null(f.g_name)
    # TO ADD: pass ahead a total number of sA that will be created by DatNet class (need it to pre-allocate self$dat.sWsA)
    # TO ADD: Current structure requires building sA twice, once for observed data and once for g_0 when g_0 unknown. This can be expensive. 
    # TO ADD: Change DatNet$new to accept vector for Odata.
    #---------------------------------------------------------------------------------
    makesWsA = function(datnetW, datnetA) {
      cbind(datnetW$dat.sVar, datnetA$dat.sVar)
    }
    gen_sWsA_dat = function(p = 1, Kmax, nodes, datnetW, datnetA, f.g_name = NULL, f.g_args = NULL)  {
      names.sWsA <- c(datnetW$names.sVar, datnetA$names.sVar)
      assert_that(is.count(p)) # self$p <- p
      nobs <- datnetW$nobs

      if (is.null(f.g_name)) {  # return observed data sW,sA if g fun is nul
        Odat.sW.sA <- makesWsA(datnetW, datnetA)
        colnames(Odat.sW.sA) <- names.sWsA
        return(data.frame(Odat.sW.sA))
      }

      dat.sWsA <- matrix(nrow = (nobs * p), ncol = (datnetW$ncols.sVar + datnetA$ncols.sVar))  # pre-allocate result matx
      colnames(dat.sWsA) <- names.sWsA

      for (i in seq(p)) {  
        Avec.df <- data.frame(Anode = f.gen.A.star(Kmax, datnetW$dat.netVar, f.g_name, f.g_args), 
                              stringsAsFactors = FALSE)
        colnames(Avec.df)[1] <- nodes$Anode
        datnetAstar <- DatNet$new(Odata = Avec.df,
                                  NetInd_k = NetInd_k, Kmax = Kmax, nodes = nodes, VarNodes = nodes$Anode)

        datnetAstar$make_sVar(names.sVar = datnetA$names.sVar) # create summary measures sA
        dat.sWsA[((i - 1) * nobs + 1):(nobs * i), ] <- makesWsA(datnetW = datnetW, datnetA = datnetAstar)
      }
      return(data.frame(dat.sWsA))
    }

    #---------------------------------------------------------------------------------
    # PARAMETERS FOR LOGISTIC ESTIMATION OF h
    #---------------------------------------------------------------------------------      
    n <- nrow(data)
    n_samp_g0gstar <- h_fit_params$n_samp_g0gstar  # replace with p adaptive to k: p <- 100*(2^k)
    family <- h_fit_params$family
    k <- h_fit_params$k
    node_l <- h_fit_params$node_l
    NetInd_k <- h_fit_params$NetInd_k
    lbound <- h_fit_params$lbound
    max_npwt <- h_fit_params$max_npwt
    logit_sep_k=h_fit_params$logit_sep_k
    h_form=h_fit_params$h_form
    f.g.star=h_fit_params$f.g.star; f.g_args=h_fit_params$f.g_args
    f.g0=h_fit_params$f.g0; f.g0_args=h_fit_params$f.g0_args
    # h_user=h_fit_params$h_user; h_user_fcn=h_fit_params$h_user_fcn; NOT IMPLEMENTED
    nFnode <- node_l$nFnode; Anode <- node_l$Anode

    # Defining the sW names that will be used for fitting h_g0/h_gstar:
    # h_form is the set of all sW (baseline summary measures) that each h_i and h^*_i depend on
    # TO DO: split into two h_forms: hform_g0 and hform_gstar? # TO BE REPLACED WITH SUMMARY MEASURES sW, sW_star, sA, sA_star (specified by the user)...
    if (is.null(h_form)) { 
      # ... NOT IMPLEMENTED .. 
      # when no regression for h/h^* is specified?
    } else {
      W_nms <- all.vars(as.formula(h_form))[-1]
    }
    if (!(nFnode%in%W_nms)) { W_nms <- c(W_nms,nFnode) }  # Always adding nFnode as a covariate

    #---------------------------------------------------------------------------------
    # BUILDING OBSERVED (netW, sW) (sW is a summary measures of netW)
    # obsdat.sW - a dataset (matrix) of n observed summary measures sW
    #---------------------------------------------------------------------------------
    # ...
    # I) Build network vectors: (W, W_netF_1, ..., W_netF_k) for each W in Wnodes by PRE-ALLOCATING netW_full:
    datnetW <- DatNet$new(Odata = data, NetInd_k = NetInd_k, Kmax = k, nodes = node_l, VarNodes = node_l$Wnodes, AddnFnode = TRUE)
    # datnetW <- DatNet$new(Odata = data, NetInd_k = NetInd_k, Kmax = k, nodes = node_l, VarNodes = node_l$Wnodes, AddnFnode = TRUE, misValRepl = TRUE)
    netW_full <- datnetW$dat.netVar
    print("datnetW$ncols.netVar: "%+%datnetW$ncols.netVar);
    # print("datnetW$names.netVar: "); print(datnetW$names.netVar)
    # ...
    # II) APPLY THE SUMMARY MEASURE FUNCTIONS / EXPRESSION TO netW_full to OBTAIN sW columns SELECT ONLY sW columns in hform_g0 and hfrom_gstar or use all? 
    sW_nms <- W_nms # change that to the actual names of summary measures in sW or entire expressions sW
    obsdat.sW <- datnetW$make_sVar(names.sVar = sW_nms)$dat.sVar
    # print("original W_nms extracted from h_form: "%+%length(W_nms)); print(W_nms);
    # print("datnetW$ncols.sVar: "%+%datnetW$ncols.sVar);
    # print("datnetW$names.sVar: "); print(datnetW$names.sVar)
    # print("creating obsdat.sW:"); print(head(obsdat.sW))
    # ...
    # III) REPLACE ALL misval values in sW with gvars$misXreplace (OR DO IT ON THE ACTUAL DATASET WHEN SUBSETTING IS PERFORMED)
    # print("replacing missing with misXreplace in obsdat.sW.");
    obsdat.sW <- datnetW$fixmiss_sVar()$dat.sVar
    #print(head(obsdat.sW))
    # ...
    # IV) (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sW:
    # print("adding determ cols to obsdat.sW (with default misval).");
    #-----------------------------------------------------------
    obsdat.sW <- datnetW$add_deterministic(Odata = data, userDETcol = "determ.g")$dat.sVar
    #print(head(obsdat.sW))
 
    #---------------------------------------------------------------------------------
    # BUILDING OBSERVED (netA, sA) (sA - summary measures of netA)
    # (actual sW_g0 and sW_gstar used for fitting h_g0 and h_gstar can be a subset of columns in sW)
    # NO NEED TO SAVE indA SEPEARTELY, AS ITS PART OF THE OBSERVED (sW,sA) and can be extracted from it when getting likelihoods
    # (netW, sW) is defined using class DatNet 
    # (netA,sA) is defined using class DatNet 
    # Final datasets for fitting h_g0/h_gstar are constructed from cbind(sW,sA), but may use netW when sampling from g0 or gstar
    #---------------------------------------------------------------------------------
    # ...
    # I) Build network vectors: (A, A_netF_1, ..., A_netF_k)
    datnetA <- DatNet$new(Odata = data, NetInd_k = NetInd_k, Kmax = k, nodes = node_l, VarNodes = node_l$Anode)
    # ...
    # II) APPLY THE SUMMARY MEASURE FUNCTIONS / EXPRESSION TO netW_full to OBTAIN sW columns
    # SELECT ONLY sW columns in hform_g0 and hfrom_gstar or use all?
    sA_nms <- netvar2(Anode, c(0:k))  # change that to the actual names of summary measures in sA or entire expressions sA
    sA_dat <- datnetA$make_sVar(names.sVar = sA_nms)$dat.sVar
    print("original sA_nms: "); print(sA_nms)
    print("datnetA$names.sVar: "); print(datnetA$names.sVar);
    print("datnetA$ncols.sVar: "%+%datnetA$ncols.sVar);
    #-----------------------------------------------------------
    # obs.sA - a dataset (matrix) of n observed summary measures sA
    obsdat.sA <- datnetA$dat.sVar # indA <- datnetA$dat.netVar
    #-----------------------------------------------------------
    # obsdat.sW.sA - a dataset (matrix) of n observed summary measures (sW,sA)
    obsdat.sW.sA <- gen_sWsA_dat(p = 1, Kmax = k, nodes = node_l, datnetW = datnetW, datnetA = datnetA)                            
    print("obsdat.sW.sA"); print(head(obsdat.sW.sA))

    #-----------------------------------------------------------
    # DEFINING SUBSETING EXPRESSIONS (FOR DETERMINISTIC / DEGENERATE sA)
    # (1 subset expr per regression P(sA[j]|sA[j-1:0], sW))
    # TO DO: Put this in a separate function (with var as arg + additional args)
    #-----------------------------------------------------------
    # (1) Capture expression as characeter string: subsetexpr <- deparse(substitute(subsetexpr))
    subsets_chr <- lapply(sA_nms, function(var) {"!misfun("%+%var%+%")"})  # subsetting by !gvars$misval on sA:
    # subsets_chr <- lapply(sA_nms, function(var) {"!misfun(nFnode)"})  # subsetting by !gvars$misval on sA:
    # subset_exprs <- lapply(sA_nms, function(var) {var%+%" != "%+%"misval"}) # compares to misval constant from gvars envir.    
    # subset_exprs <- lapply(netvar2("determ.g_Friend", c(0:k)), function(var) {var%+%" != "%+%"misval"}) # based on the variable of gvars$misval (requires passing gvars envir for eval)
    # subset_exprs <- lapply(netvar2("determ.g_true", c(0:k)), function(var) {var%+%" != "%+%TRUE}) # based on existing logical determ_g columns (TRUE = degenerate/determ):
    # (2) Parse the characteer expression into call (make subset expressions into a list of calls (one call per sA[j] in sA))
    # (3) Substitute the actual var names in the data for generic node names (nFnode, Wnodes, Anode):
    substitute_list <- lapply(node_l, as.name)
    subsets_expr <- lapply(subsets_chr, function(subset_chr) {      
                                            subset_expr <- try(parse(text=subset_chr)[[1]]) # parses chr into a call
                                            if(inherits(subset_expr, "try-error")) stop("can't parse the subset formula", call.=FALSE)
                                            eval(substitute(substitute(e, env = substitute_list), list(e = subset_expr)))
                                          })
    # print("subsets_expr: "); print(subsets_expr)

    ##########################################
    # Summary class params:
    ##########################################
    # sA_class <- c("binary", "contin", rep_len("binary", 5))
    sA_class <- rep_len("binary", length(sA_nms))

    summary_params <- list(sA_class = sA_class, sA_nms = sA_nms, sW_nms = W_nms, subset = subsets_expr)
    # summary_params <- list(Kmax = k, nodes = node_l, sA_class = sA_class, sA_nms = sA_nms, sW_nms = W_nms, subset = subsets_expr)
    ##########################################
    message("fitting h under g_0...")
    ##########################################
    p_h0 <- ifelse(is.null(f.g0), 1, n_samp_g0gstar)
    fit.g0_dat <- gen_sWsA_dat(p = p_h0, Kmax = k, nodes = node_l, datnetW = datnetW, datnetA = datnetA, 
                                f.g_name = f.g0, f.g_args = f.g0_args)
    # print("fit.g0_dat: "); print(head(fit.g0_dat))
    # above fun gen_sWsA_dat will becomes part of DatSummaries class:
    # DatSummaries$new(Odata = data, Kmax = k, nodes = node_l, NetInd_k = NetInd_k, f.g_name = f.g0, f.g_args = f.g_args)

    summeas.g0 <- do.call(SummariesModel$new, summary_params)
    # print("summeas.g0$regs_list: "); print(summeas.g0$regs_list)
    summeas.g0$fit(data = fit.g0_dat)
    summeas.g0$predict(newdata = obsdat.sW.sA) # DOESN'T HAVE TO BE CALLED IF (is.null(f.g0)), since PREDICATIONS ARE ALREADY SAVED for obsdat.sW.sA:
    summeas.g0$predictAeqa(obsdat.sA = obsdat.sA)
    h_vec.g0.new <- summeas.g0$getcumprodAeqa()

    ##########################################
    message("fitting h under g_star...")
    ##########################################
    fit.gstar_dat <- gen_sWsA_dat(p = n_samp_g0gstar, Kmax = k, nodes = node_l, datnetW = datnetW, datnetA = datnetA, 
                                    f.g_name = f.g.star, f.g_args = f.g_args)    
    # print("fit.gstar_dat: "); print(head(fit.gstar_dat)); 
    
    summeas.gstar <- do.call(SummariesModel$new, summary_params)
    summeas.gstar$fit(data = fit.gstar_dat)
    summeas.gstar$predict(newdata = obsdat.sW.sA)
    summeas.gstar$predictAeqa(obsdat.sA = obsdat.sA)
    h_vec.gstar.new <- summeas.gstar$getcumprodAeqa()

    ###########################################
    # alternative to above using a function call instead
    ###########################################
    # fit.and.predict.h.new <- function(Kmax, sA_nms, sW_nms, hfitdat, newdata, obsdat.sA) { # NOT USED FOR NOW:
    #   summeas.g <- SummariesModel$new(Kmax = Kmax, sA_nms = sA_nms, sW_nms = W_nms)
    #   summeas.g$fit(data = hfitdat)
    #   if (!missing(newdata)) { # don't need to predict again if only need predictions for fit data
    #     summeas.g$predict(newdata = newdata)
    #   }
    #   summeas.g$predictAeqa(obsdat.sA = obsdat.sA)
    #   h_vec.g <- summeas.g$getcumprodAeqa()
    #   return(list(h_vec.g = h_vec.g, summeas.g = summeas.g))
    # }
    # h_g0.new <- fit.and.predict.h.new(Kmax = k, sA_nms = sA_nms, sW_nms = W_nms, hfitdat = fit.g0_dat, obsdat.sA = obsdat.sA)  # new method based on R6 classes:
    # summeas.g0 <- h_g0.new$summeas.g
    # h_vec.g0.new <- h_g0.new$h_vec.g
    # h_gstar.new <- fit.and.predict.h.new(Kmax = k, sA_nms = sA_nms, sW_nms = W_nms, hfitdat = fit.gstar_dat, newdata = obsdat.sW.sA, obsdat.sA = obsdat.sA)  # new method based on R6 classes:
    # summeas.gstar <- h_gstar.new$summeas.g
    # h_vec.gstar.new <- h_gstar.new$h_vec.g


    ###########################################
    # 3) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
    ##########################################
    h_tilde.new <- h_vec.gstar.new / h_vec.g0.new
    h_tilde.new[is.nan(h_tilde.new)] <- 0     # 0/0 detection
    h_tilde.new <- bound(h_tilde.new, c(0,1/lbound))

    df_h_bar_vals <- data.frame(cY.ID = .f.mkstrNet(obsdat.sW.sA), 
                                h.star_c = h_vec.gstar.new,
                                h_c = h_vec.g0.new,
                                h = h_tilde.new
                                )
    print("predicted h for obs. data:"); print(head(df_h_bar_vals, 20))

    fit_h_reg_obj <- list(k=k,
                          m.gAi_vec_g = summeas.g0,
                          m.gAi_vec_gstar = summeas.gstar,
                          lbound=lbound,
                          # determ_cols_Friend=determ_cols_Friend, # determ_cols_fitted=determ_cols, 
                          cY_mtx_fitted=obsdat.sW.sA
                          )

    return(list(df_h_bar_vals=df_h_bar_vals, fit_h_reg_obj=fit_h_reg_obj))
}

# #---------------------------------------------------------------------------------
# # Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
# #---------------------------------------------------------------------------------
# .f_get_all_ests <- function(data, data_net, est_obj, MCeval_hstar) {
#   n <- nrow(data)
#   node_l <- est_obj$node_l
#   nFnode <- node_l$nFnode
#   Anode <- node_l$Anode
# 	Ynode <- node_l$Ynode
# 	Y <- data[, Ynode]
# 	determ.Q <- data[, "determ.Q"]
#   determ.g <- data[, "determ.g"]

#   #************************************************
#   # IPTW_h estimator (based on h^*/h_N clever covariate):
#   #************************************************
#   fit.hbars_t <- system.time(h_bars <- fit.hbars(data=data, h_fit_params=est_obj)) # fit the clever covariat
#   # print("time to fit h_bars"); print(fit.hbars_t)
#   df_h_bar_vals <- h_bars$df_h_bar_vals
#   fit_h_reg_obj <- h_bars$fit_h_reg_obj
#   h_wts <- df_h_bar_vals$h
#   Y_h_wts <- Y
#   Y_h_wts[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]
#   h_iptw <- mean(Y_h_wts)  # IPTW estimator based on h - clever covariate
#   # print("IPW Est (h)"); print(mean(Y_h_wts))

#   #************************************************
#   # IPTW_g estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
#   #************************************************
# 	# 02/16/13: IPTW estimator (Y_i * prod_{j \in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
# 	g_wts <- iptw_est(k=est_obj$k, data=data, node_l=node_l, m.gN=est_obj$m.g0N, f.g.star=est_obj$f.g.star, f.g_args=est_obj$f.g_args, family=est_obj$family, 
#                       NetInd_k=est_obj$NetInd_k, lbound=est_obj$lbound, max_npwt=est_obj$max_npwt, f.g0=est_obj$f.g0, f.g0_args=est_obj$f.g0_args)
#   Y_g_wts <- Y
#   Y_g_wts[!determ.Q] <- Y[!determ.Q] * g_wts[!determ.Q]
#   g_iptw <- mean(Y_g_wts)  # IPTW estimator based on full g factorization (prod(g))

# 	#-------------------------------------------
#   return(list( h_iptw = h_iptw,
#                g_iptw = g_iptw,
#                h_wts=h_wts, g_wts=g_wts))
# }

# #---------------------------------------------------------------------------------
# # MAIN TMLE ESTIMATOR FUNCTION
# #---------------------------------------------------------------------------------
# tmlenet <- function(data, Anode, Wnodes, iidW_flag=FALSE, Ynode, nFnode, 
#                     k_max, IDnode=NULL, NETIDnode,
#                     Qform=NULL, QDETnode=NULL,
#                     Q.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
#                     gform=NULL, gDETnode=NULL, gbound=0.005, max_npwt=50, 
#                     g.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
#                     f.g1.star, f.g1_args, f.g2.star=NULL, f.g2_args=NULL, 
#                     h_f.g0=NULL, h_f.g0_args=NULL, h_user_fcn = NULL, h_form = NULL,
#                     h_logit_sep_k = FALSE, W_indep=FALSE,  family = "binomial", 
#                     alpha  = 0.05, verbose=FALSE, n_MCsims=ceiling(sqrt(nrow(data))), n_samp_g0gstar=20, 
#                     alphaTMLE_B_only = TRUE) {

#   # if (is.na(h_form)) h_form <- NULL
#   print("Running tmlenet with "); 
#   print("Qform"); print(Qform); 
#   print("gform"); print(gform); 
#   print("hform"); print(h_form);

#   # DEBUGGING:
#   # memory & speed profiling
#   # Rprof("tmle_run_memuse.out", memory.profiling=TRUE)
#   #------------------------------------------------------------------------------
#   # Which esimators to evaluate
#   #------------------------------------------------------------------------------
#   MCeval_hstar <- FALSE   # if FALSE, do not evaluate MC wrt density h_star
#   #------------------------------------------------------------------------------
#   # MONTE-CARLO SIMULATION PARAMETERS
#   #------------------------------------------------------------------------------
#   n_MCsims <- as.integer(n_MCsims)  # number of times to sample for MC sim  
#   max.err_est <- 0.1    # maximum percent error for MCS estimators
#   #------------------------------------------------------------------------------
#   # ESTIMATION OF h & h^* PARAMETERS
#   #------------------------------------------------------------------------------
#   n_samp_g0gstar <- n_samp_g0gstar # number of times to sample from g^* or g_0 (if known) for estimation of \bar{h}^* and \bar{h}_0
# 	logit_sep_k = h_logit_sep_k # Fit each k value separately (T) or all at once (F) (now passed as arguments to tmlenet)
#   f.g0 <- h_f.g0
#   f.g0_args <- h_f.g0_args
# 	d <- data
# 	n <- nrow(d)
# 	k <- k_max
#   h_user <- !(is.null(h_user_fcn))

#   # **********
#   #Q: 02/15/14: IS IT A GOOD IDEA TO HAVE THE SAME (OR ANY) UPPER BOUND ON g_star? Doesn't make any sense...
#   # **********
# 	if(length(gbound)==1) gbound <- c(gbound, 1-gbound)
#  	#------------------------------------------------------------------------------
# 	# Netwk ids strings to list of friend indices (NetVec) and matrix of friend indices (NetInd_k)
#   .f_getNetIndices <- function (d) {
#     .f.mkvecNetID <- function(Net_str) lapply(Net_str, function(Net_str_i) unlist(strsplit(Net_str_i, ' ', fixed=TRUE)))
#       # Netwk ids strings to arr of ID vectors (filled up with trailing NA's)
#     .f.mkarrNetID <- function(Net_str) t(sapply(Net_str, function(Net_str_i) {
#                       netwk <- unlist(strsplit(Net_str_i, ' ',
#                       fixed=TRUE))
#                       return(c(netwk, rep_len(NA,k-length(netwk))))
#                       } ))
#     NetIDVec <- .f.mkvecNetID(d[,NETIDnode])  # get list of network ID vectors from network strings for each i
#     NetVec <- lapply(NetIDVec, function(NetID) as.numeric(sapply(NetID, function(x) which(x==as.vector(d[,IDnode]))))) # convert into list of network row #s
#     NetID_k <- .f.mkarrNetID(d[,NETIDnode]) # get array of network IDs
#     NetInd_k <- apply(NetID_k, 2, function(k_ID) {  # make NetID_k into array of network indices (row #s)
#                   sapply(as.vector(k_ID), function(x) {
#                     if (is.na(x)) {
#                       NA
#                     } else { 
#                       which(x==as.vector(d[,IDnode]))
#                     }
#                    })
#                 })
#     return(list(NetVec=NetVec, NetInd_k=NetInd_k))
#   }
#   NetInd_l <- .f_getNetIndices(d)
#   NetVec <- NetInd_l$NetVec
#   NetInd_k  <- NetInd_l$NetInd_k
#   #------------------------------------------------------------------------------
#   # List of node names and networks of W's and A's (netA and netW)
#   #------------------------------------------------------------------------------
#   node_l <- list(IDnode=IDnode, Anode=Anode, Wnodes=Wnodes, Ynode=Ynode, nFnode=nFnode, NETIDnode=NETIDnode)
#   if (is.null(gDETnode)) determ.g <- rep_len(FALSE,n) else 
#       determ.g <- (d[, gDETnode] == 1)
#   if (is.null(QDETnode)) determ.Q <- rep_len(FALSE,n) else 
#       determ.Q <- (d[, QDETnode] == 1)

# 	netW <- NULL
# 	for (Wnode in node_l$Wnodes) {
# 	 	netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, d[,Wnode], Wnode)))
# 	}
# 	netA <- data.frame(.f.allCovars(k, NetInd_k, d[,node_l$Anode], node_l$Anode))
# 	df_AllWs <- netW	 
# 	net_d <- cbind(ID=d[, node_l$IDnode], netW, netA, subset(d, select=c(node_l$nFnode,node_l$Ynode)))
# 	# print("net_d"); print(head(net_d, 10))

# 	#-------------------------------------------					
# 	# Fit initial model for Q(Y|A,W) under oberved data (m.Q.init) - move this to .f_get_all_ests() to avoid confusion
# 	#-------------------------------------------  
# 	d_sel <- data.frame(subset(d, select=unlist(node_l)), determ.g=determ.g, determ.Q=determ.Q, QY.init=QY.init, iidQY.init=iidQY.init)
#   node_l <- c(node_l, gform=gform, Qform=Qform, iidW_flag=iidW_flag)
#   # print("m.Q.init"); print(summary(m.Q.init)); print("QY.init fit"); print(QY.init);
# 	#-------------------------------------------							
# 	# Fit the model for g_0(A,W) - move this to .f_get_all_ests() to avoid confusion
# 	#-------------------------------------------
# 	m.g0N <- .f.est(net_d[!determ.g,], gform, family=family) # Set A=0 when determ.g==1
#   #-------------------------------------------              
#   # Create an object with model estimates, data & network information that is passed on to estimation procedure
#   #-------------------------------------------
#   # 1) define parameters for MC estimation of the substitution estimators
#   # 2) define parameters for estimation of the efficient weights h(A^s|W^s)
#   est_obj <- list(
#     k=k, node_l=node_l, NetInd_k=NetInd_k,
#     lbound=gbound[1], 
#     family = family,
#     m.g0N=m.g0N,
#     f.g0=f.g0, f.g0_args=f.g0_args,
#     logit_sep_k=logit_sep_k, 
#     h_user=h_user, 
#     h_user_fcn=h_user_fcn, 
#     h_form=h_form, 
#     n_samp_g0gstar=n_samp_g0gstar)
	
#   #-------------------------------------------							
# 	# Run TMLE for each g.star and/or ATE
# 	#-------------------------------------------
#   est_obj_g1 <- append(est_obj, list(f.g.star=f.g1.star, f.g_args=f.g1_args))
# 	tmle_g1_out <- .f_get_all_ests(data=d_sel, data_net=net_d, est_obj=est_obj_g1, MCeval_hstar=MCeval_hstar)
# 	if (!is.null(f.g2.star)) {
#     est_obj_g2 <- append(est_obj, list(f.g.star=f.g2.star, f.g_args=f.g2_args))
# 		tmle_g2_out <- .f_get_all_ests(data=d_sel, data_net=net_d, est_obj=est_obj_g2, MCeval_hstar=MCeval_hstar)
# 	}
# 	else {
# 		tmle_g2_out <- NULL
# 	}

#   # create output object with param ests of EY_gstar, vars and CIs for given gstar (or ATE if two tmle obj are passed)
#   .f_make_EYg_obj <- function(tmle_g_out, tmle_g2_out=NULL) {
#     # get estimates of as. var and CIs:
#     if (is.null(tmle_g2_out)) {
#       tmle_g2_out <- list()
#       tmle_g2_out$QY.star <- tmle_g2_out$fWi_init_A <- tmle_g2_out$fWi_init_B <- tmle_g2_out$fWi_star_A <- tmle_g2_out$fWi_star_B <- tmle_g2_out$fWi_init_tmleiptw <- tmle_g2_out$h_iptw <- tmle_g2_out$iptw_reg <- 0
#       tmle_g2_out$tmle_A <- tmle_g2_out$tmle_B <- tmle_g2_out$iid.tmle_B <- tmle_g2_out$tmle_iptw <- tmle_g2_out$iptw_h <- tmle_g2_out$iptw <- tmle_g2_out$iid.iptw <- tmle_g2_out$mle <- 0
#       tmle_g2_out$noMC_tmle_A <- tmle_g2_out$noMC_tmle_B <- tmle_g2_out$noMC_mle <- 0
#     }

#     tmle_A <- tmle_g_out$tmle_A - tmle_g2_out$tmle_A
#     noMC_tmle_A <- tmle_g_out$noMC_tmle_A - tmle_g2_out$noMC_tmle_A

#     tmle_B <- tmle_g_out$tmle_B - tmle_g2_out$tmle_B
#     iid.tmle_B <- tmle_g_out$iid.tmle_B - tmle_g2_out$iid.tmle_B
#     noMC_tmle_B <- tmle_g_out$noMC_tmle_B - tmle_g2_out$noMC_tmle_B

#     tmle_iptw <- tmle_g_out$tmle_iptw - tmle_g2_out$tmle_iptw
#     iptw_h <- tmle_g_out$iptw_h - tmle_g2_out$iptw_h
#     iptw <- tmle_g_out$iptw - tmle_g2_out$iptw
#     iid.iptw <- tmle_g_out$iid.iptw - tmle_g2_out$iid.iptw
    
#     mle = tmle_g_out$mle - tmle_g2_out$mle
#     noMC_mle = tmle_g_out$noMC_mle - tmle_g2_out$noMC_mle
 
#     EY_g.star <- list(  tmle_A = as.vector(tmle_A),
#                         noMC_tmle_A = as.vector(noMC_tmle_A),
#                         tmle_B = as.vector(tmle_B),
#                         iid.tmle_B = as.vector(iid.tmle_B),
#                         noMC_tmle_B = as.vector(noMC_tmle_B),
#                         tmle_iptw = as.vector(tmle_iptw),
#                         iptw_h = as.vector(iptw_h),
#                         iptw = as.vector(iptw),
#                         iid.iptw = as.vector(iid.iptw),
#                         mle = as.vector(mle),
#                         noMC_mle = as.vector(noMC_mle)
#                         )
#     return(EY_g.star)
#   }

#   EY_g1.star <- .f_make_EYg_obj(tmle_g_out=tmle_g1_out)

# 	if (!is.null(f.g2.star)) {	
#     EY_g2.star <- .f_make_EYg_obj(tmle_g_out=tmle_g2_out)
#     ATE <- .f_make_EYg_obj(tmle_g_out=tmle_g1_out, tmle_g2_out=tmle_g2_out)
# 	} else {
# 		EY_g2.star <- NULL
# 		ATE <- NULL
# 	}
#   # print("EY_g1.star"); print(EY_g1.star)
#   # print("EY_g2.star"); print(EY_g2.star)
# 	estimates <- list(EY_g1.star=EY_g1.star, EY_g2.star=EY_g2.star, ATE=ATE)
# 	tmlenet <- list(estimates=estimates)
# 	class(tmlenet) <- "tmlenet"	

#   out_sim <- rbind(EY_g1.star$tmle_A, EY_g1.star$tmle_B, EY_g1.star$iid.tmle_B, EY_g1.star$tmle_iptw, EY_g1.star$iptw_h, EY_g1.star$iptw, EY_g1.star$iid.iptw, EY_g1.star$mle)
#   rownames(out_sim) <- c("tmle_A", "tmle_B", "iid.tmle_B", "tmle_iptw", "iptw_h", "iptw", "iid.iptw", "mle")
#   print("Estimates w/ MC eval:"); print(out_sim)

#   noMC_out_sim <- rbind(EY_g1.star$noMC_tmle_A, EY_g1.star$noMC_tmle_B, EY_g1.star$noMC_mle)
#   rownames(noMC_out_sim) <- c("noMC_tmle_A","noMC_tmle_B","noMC_mle")
#   # print("Estimates w/out MC eval:"); print(noMC_out_sim)

#   # DEBUGGING:
#   # Rprof(NULL)
# 	return(tmlenet)
# }
# #----------------------------------------------------------------------------------	  
# # #   browser {base}
# # R Documentation
# # Environment Browser
# # Description
# # Interrupt the execution of an expression and allow the inspection of the environment where browser was called from.
# #----------------------------------------------------------------------------------	  