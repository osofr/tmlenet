# #todo 61 (MCeval) +0: Consider adding model update fits to mcEvalPsi, as separate methods or as part of each TMLE

#---------------------------------------------------------------------------------
# G-Comp & TMLEs: Use Monte-Carlo to estimate psi under stochastic g^* 
#---------------------------------------------------------------------------------
get.MCS_ests <- function(datNetObs,  DatNet.gstar, MC_fit_params, m.h.fit) {
  onlyTMLE_B <- MC_fit_params$onlyTMLE_B
  # nQ.MCsims <- MC_fit_params$nQ.MCsims
  # max_npwt <- MC_fit_params$max_npwt
  # max.err_eps <- MC_fit_params$max.err_eps
  # f.gstar <- MC_fit_params$f.gstar
  # f.g0 <- MC_fit_params$f.g0
  m.Q.init <- MC_fit_params$m.Q.init
  m.Q.starA <- MC_fit_params$m.Q.starA
  m.Q.starB <- MC_fit_params$m.Q.starB
  # m.Q.star_giptw <- MC_fit_params$m.Q.star_giptw
  # m.g0N  <- MC_fit_params$m.g0N

  message("================================================================")
  message("running Monte Carlo evaluation of the substitution estimators...")
  message("================================================================")
  evaluator <- mcEvalPsi$new(datNetObs = datNetObs, DatNet.gstar = DatNet.gstar, onlyTMLE_B = onlyTMLE_B)
  nOdata <- evaluator$nOdata

  genMC.reps <- function(nrep)  {
    # resamp_d <- .f.gen.sample(NetInd_k, n) # Get a random sample of all A and W (DISABLED)
    GCOMP <- evaluator$get.gcomp(m.Q.init) # QY.init (G-Comp estimator) - est probY based on model for Q_Y

    TMLE_A <- 0 # TMLE A (adjusted by coefficient epsilon on h_bar ratio)
    TMLE_gIPTW <- 0 # gIPTW  TMLE
    if (!onlyTMLE_B) {
      TMLE_A <- evaluator$get.tmleA(m.Q.starA = m.Q.starA, m.h.fit = m.h.fit) # TMLE A (adjusted by coefficient epsilon on h_bar ratio)
      TMLE_gIPTW <- 0 # gIPTW  TMLE
    }

    TMLE_B <- evaluator$get.tmleB(m.Q.starB = m.Q.starB) # TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
    fiWs_list <- evaluator$get.fiW(m.Q.starA = m.Q.starA, m.Q.starB = m.Q.starB) # Get fi_W - hold W fixed to observed values

    # Put all estimators together and add names (defined in out_nms outside of this function):
    mean_psis_all <- c(tmle_A = mean(TMLE_A), tmle_B = mean(TMLE_B), tmle_g_iptw = mean(TMLE_gIPTW), mle = mean(GCOMP),
                      fiWs_list$fiW_Qinit, 
                      fiWs_list$fiW_QstarA, 
                      fiWs_list$fiW_QstarB)

    # Naming estimators (same names throughout the entire package):
    names(mean_psis_all) <- c("tmle_A", "tmle_B", "tmle_g_iptw", "mle", 
                            paste("fWi_init_", c(1:nOdata), sep = ""),
                            paste("fWi_star_A_", c(1:nOdata), sep = ""),
                            paste("fWi_star_B_", c(1:nOdata), sep = ""))

    return(mean_psis_all)
  }

  psi_est_mean <- genMC.reps(1)

  # ********************************************************
  # HAVE NOT YET IMPLEMENTED DYNAMIC MC CONVERGENCE TESTING IN THIS VS:
  # ********************************************************
  #---------------------------------------------------------------------------------
  # Main body of a fcn get.MCS_ests(): MC evalution of the estimators
  #---------------------------------------------------------------------------------
  # # Allow this part to loop, until desired MCS prob_epsilon for all estimators is reached:
  # nrepeat <- 1
  # psis_reps <- NULL
  # ests_reps <- NULL
  # repeat {
  #   ests_reps <- rbind(ests_reps,
  #                       t(sapply(seq(nQ.MCsims) , genMC.reps))
  #                       )
  #   psi_est_mean <- apply(ests_reps, 2, mean, na.rm = T)
  #   psi_est_var <- apply(ests_reps, 2, var, na.rm = T)
  #   psi_percerr <- 2 * abs(psi_est_mean * max.err_eps) # estimate the maximum allowed epsilon for each estimator, based pre-defined % error:  
  #   prob_percerr <- psi_est_var / ((nQ.MCsims * nrepeat) * (psi_percerr)^2)
  #   prob_percerr[psi_est_var < 0.0001] <- 0.0001
  #   fin_ests_sel <- c(1:3) # final vec of estimators for which error is measured
  #   if ( (all(prob_percerr[fin_ests_sel] < 0.05)) | (nrepeat >= 100)) {
  #     break
  #   }
  #   nrepeat <- nrepeat + 1
  # }   
  # # print("nrepeat"); print(nrepeat)

  return(psi_est_mean)
}


# Class for Monte-Carlo evaluation of various substitution estimators under user-specified stochastic intervention gstar
  # For given data, take Q[Y|cY]=m.Q.init and calcualte est. of psi under gstar using Monte-Carlo integration:
  # * W_i can be iid or not (W's are never resampled);
  # * Draw from the distributions of W and \bar{g}^*
  # * Recalculate Y^c under \bar{g}^*;
  # * Repeat nrep times and average.
#' @import data.table

## ---------------------------------------------------------------------
#' @title Class for Monte-Carlo evaluation of various substitution estimators under user-specified stochastic intervention gstar.
#' @docType class
#' @format An R6 class object.
#' @name mcEvalPsi
#' @details Following fields are created during initialization
#' \itemize{
#' \item{nodes} ...
#' \item{subset_regs} ...
#' \item{sA_nms} ...
#' \item{sW_nms} ...
#' \item{Kmax} ...
#' }
#' Class for evaluating and storing arbitrary summary measures sVar.
#' @importFrom assertthat assert_that is.count is.flag
# @export
mcEvalPsi <- R6Class(classname = "mcEvalPsi",
  portable = TRUE,
  class = TRUE,
  public = list(
    Kmax = integer(),          # max n of Friends in the network
    nodes = list(),            # names of the nodes in the data (Anode, Wnodes, nFnode, etc..)
    netind_cl = NULL,          # class NetIndClass object holding $NetInd_k network matrix

    datNetObs = NULL,
    DatNet.gstar = NULL,

    onlyTMLE_B = NULL,
    m.Q.init = NULL,
    m.Q.starA = NULL,
    m.Q.starB = NULL,
    QY.init = NULL,
    QY.starA = NULL,
    QY.starB = NULL,

    nOdata = NA_integer_,      # no. of samples in the OBSERVED (original) data
    p = NA_integer_,           # no. of times n-size sA were resampled under gstar

    initialize = function(datNetObs, DatNet.gstar, onlyTMLE_B, ...) {
      self$datNetObs <- datNetObs
      self$nOdata <- datNetObs$nobs
      self$Kmax <- datNetObs$Kmax
      self$nodes <- datNetObs$nodes
      self$netind_cl <- datNetObs$netind_cl
      self$DatNet.gstar <- DatNet.gstar
      self$p <- DatNet.gstar$p

      # Additional params:
      self$onlyTMLE_B <- onlyTMLE_B
      # nQ.MCsims <- MC_fit_params$nQ.MCsims
      # max_npwt <- MC_fit_params$max_npwt
      # max.err_eps <- MC_fit_params$max.err_eps

      # f.gstar <- MC_fit_params$f.gstar
      # f.g0 <- MC_fit_params$f.g0
      # m.g0N  <- MC_fit_params$m.g0N
      invisible(self)
    },

    # MLE - Predict E[Yg_star] (QY.init) for each i, based on the initial model fit for E[Y|C^Y] (m.Q.init)
    # output is a vector of length n*p
    get.gcomp = function(m.Q.init) {
      # QY <- predict(m.Q.init, newdata = samp_data, type="response")
      # QY[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect if W's are resampled
      self$m.Q.init <- m.Q.init
      datNetObs <- self$datNetObs
      DatNet.gstar <- self$DatNet.gstar
      QY.init <- m.Q.init$predict(newdata = DatNet.gstar)$getprobA1
      QY.init[datNetObs$det.Y] <- datNetObs$noNA.Ynodevals[datNetObs$det.Y]
      self$QY.init <- QY.init
      invisible(QY.init)
    },

    # TMLE A - Update QY.init based on the est. coefficient for the clever covariate h_g0/h_gstar in univar. model fluct (m.Q.star_h_A)
    # output is a vector of length n*p
    get.tmleA = function(m.Q.starA, m.h.fit) {
      if (is.null(self$QY.init)) stop("call mcEvalPsi$get.gcomp first") # QY.init <- self$get.gcomp(self$m.Q.init)
      predh.res <- pred.hbars(newdatnet = self$DatNet.gstar, m.h.fit = m.h.fit)
      h_iptw <- predh.res$dat_hest$h_gstar_gN
      if (!is.na(coef(m.Q.starA))) {
        self$m.Q.starA <- m.Q.starA
        off <- qlogis(self$QY.init)
        self$QY.starA <- plogis(off + coef(m.Q.starA)*h_iptw)
      } else {
        self$QY.starA <- self$QY.init
      }
      # QY.starA[determ.Q] <- samp_data[determ.Q, Ynode] # will be incorrect if W's are resampled
      invisible(self$QY.starA)
    },

    # TMLE B - Update QY.init based on the est. intercept of the model fluct (m.Q.star_h_B)
    # output is a vector of length n*p
    get.tmleB = function(m.Q.starB) {
      if (is.null(self$QY.init)) stop("call mcEvalPsi$get.gcomp first") # QY.init <- self$get.gcomp(self$m.Q.init)
      if (!is.na(coef(m.Q.starB))) {
        self$m.Q.starB <- m.Q.starB
        off <- qlogis(self$QY.init)
        self$QY.starB <- plogis(off + coef(m.Q.starB))
      } else {
        self$QY.starB <- self$QY.init
      }
      # QY.star_B[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect when W's are resampled
      invisible(self$QY.starB)
    },

    # get an estimate of fiW (hold ALL W's fixed at once) - a component of TMLE Var
    # Creates a vector of size n*p, where each of n obs is then averaged p times.
    get.fiW = function(m.Q.starA, m.Q.starB) {
      if (is.null(self$QY.init)) stop("call mcEvalPsi$get.gcomp first") # QY.init <- self$get.gcomp(self$m.Q.init)

      # *******fi_W based on Q,N.init model ******
      ID <- rep.int(c(1 : self$nOdata), self$p)
      # taking average over p samples for each of n obs
      fiW_Qinit <- data.table(ID = ID, fiW = self$QY.init)
      fiW_Qinit.mean <- fiW_Qinit[, lapply(.SD, mean, na.rm=TRUE), by="ID", .SDcols=c("fiW") ][["fiW"]]
      fiW_Qinit.var <- fiW_Qinit[, lapply(.SD, var, na.rm=TRUE), by="ID", .SDcols=c("fiW") ][["fiW"]]

      # *******fi_W based on Q.N.star models (A & B) ******
      fiW_QstarA <- rep_len(0, self$nOdata)
      fiW_QstarB <- rep_len(0, self$nOdata)

      # if (!self$onlyTMLE_B) {
      #   fiW_QstarA <- self$get.tmleA(m.Q.starA) # or: fiW_QstarA <- self$QY.starA
      # } else {
      #   fiW_QstarA <- rep_len(0, self$nOdata)
      # }
      # fiW_QstarB <- self$get.tmleB(m.Q.starB) # or: fiW_QstarB <- self$QY.starB

      # #todo 60 (MCeval, get.fiW) +0: Consider doing ONE type of fi_W at a time to simplify (based on one m.Q.init, m.Q.starA, m.Q.starB that is passed as arg)
      # return(list(fiW_Qinit = fiW_Qinit.mean, fiW_Qinit.var = fiW_Qinit.var, fiW_QstarA = fiW_QstarA, fiW_QstarB = fiW_QstarB))
      return(list(fiW_Qinit = fiW_Qinit.mean, fiW_QstarA = fiW_QstarA, fiW_QstarB = fiW_QstarB))
    }

    # NOT IMPLEMENTED
    # # IPTW NETWORK TMLE
    # get.TMLE.gIPTW <- function(QY.init) {
    #   g_iptw <- iptw_est(k=k, data=samp_data, node_l=node_l, m.gN=m.g0N, f.gstar=f.gstar, f.g_args=f.g_args,
    #                       family=family, NetInd_k=NetInd_k, lbound=lbound, max_npwt=max_npwt, f.g0=f.g0, f.g0_args=f.g0_args)
    #   determ.Q <- samp_data$determ.Q
    #   if (!is.na(coef(m.Q.star_iptw))) {
    #     off <- qlogis(QY.init)
    #     QY.star <- plogis(off + coef(m.Q.star_iptw)*g_iptw)
    #     #************************************************
    #     # QY.star[determ.Q] <- 1
    #     QY.star[determ.Q] <- samp_data[determ.Q, Ynode]
    #     #************************************************
    #     return(QY.star)
    #   }
    # },
  ),

  active = list(
    dat.sVar = function() { self$mat.sVar },
    emptydat.sVar = function() { self$mat.sVar <- NULL }       # wipe out mat.sVar
  ),
  private = list(
    placeholder = list()
  )
)