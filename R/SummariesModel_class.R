#' @importFrom assertthat assert_that

## **********************************************************************
## DESCRIPTION OF THE FITTING ALGORITHM(s) FOR CONTINUOUS SUMMARY MEASURE sA (TO ADD TO THE PACKAGE MANUAL)
## LEFT TO DO: 
# *) implement the hazard weighting by bin width for pooled cont sA fitting (in BinDat$logispredict.long)
## **********************************************************************
# *) Fit one density for P_{g_0}(sA | sW) implied by (A,W)~g_0(A|W)Q_0(W) 
#    and another density for P_{g_star}(sA^* | sW) that is implied by (A^*,W)~g_star(A|W)Q_0(W).
# *) Same algorithm is used for estimation of P_{g_0}(sA | sW) or P_{g_star}(sA^* | sW).
# *) These two densities form the basis of the IPTW estimator,
#    which is evaluated at the observed N data points o_i=(y_i, sa_i, sw_i), i=1,...,N.
#    The IPTW is then given by \sum_i={1,...,N} {Y_i * P_{g_star}(sA^*=sa_i | sW=sw_i) / P_{g_0}(sA=sa_i | sW=sw_i)}
# *) For simplicity, we now suppose sA is univariate and we first describe an algorithm for fitting P_{g_0}(sA | sW):
  # 1) Generate a dataset of N observed continuous summary measures (sa_i:i=1,...,N) from observed ((a_i,w_i):i=1,...,N). Let sa\in{sa_i:i=1,...,M}.
  # 2) Divide the range of sA values into intervals S=(i_1,...,i_M,i_{M+1}) so that any observed data point sa_i belongs to one interval in S, namely, 
  #    for each possible value sa of sA there is k\in{1,...,M}, such that, i_k < sa <= i_{k+1}.
  #    Let the mapping B(sa)\in{1,...,M} denote a unique interval in S for sa, such that, i_{B(sa)} < sa <= i_{B(sa)+1}.
  #    Let bw_{B(sa)}:=i_{B(sa)+1}-i_{B(sa)} be the length of the interval (bandwidth) (i_{B(sa)},i_{B(sa)+1}).
  #    Also define the binary indicators b_1,...,b_M, where b_j:=I(B(sa)=j), for all j <= B(sa) and b_j:=NA for all j>B(sa).
  #    That is we set b_j to missing ones the indicator I(B(sa)=j) jumps from 0 to 1.
  #    Now let sA denote the random variable for the observed summary measure for one unit
  #    and denote by (B_1,...,B_M) the corresponding random indicators for sA defined as B_j := I(B(sA) = j) 
  #    for all j <= B(sA) and B_j:=NA for all j>B(sA).
  # 3) For each j=1,...,M, fit the logistic regression model for the conditional probability P(B_j = 1 | B_{j-1}=0, sW), i.e., 
  #    at each j this is defined as the conditional probability of B_j jumping from 0 to 1 at bin j, given that B_{j-1}=0 and 
  #    each of these logistic regression models is fit only among the observations that are still at risk of having B_j=1 with B_{j-1}=0.
  # 4) Normalize the above conditional probability of B_j jumping from 0 to 1 by its corresponding interval length (bandwidth) bw_j to 
  #    obtain the discrete conditional hazards h_j(sW):=P(B_j = 1 | (B_{j-1}=0, sW) / bw_j, for each j.
  #    For the summary measure sA, the above conditional hazard h_j(sW) is equal to P(sA \in (i_j,i_{j+1}) | sA>=i_j, sW), 
  #    i.e., this is the probability that sA falls in the interval (i_j,i_{j+1}), conditional on sW and conditional on the fact that
  #    sA does not belong to any intervals before j.
  # 4) Finally, for any given data-point (sa,sw), evaluate the discretized conditional density for P(sA=sa|sW=sw) by first 
  #    evaluating the interval number k=B(sa)\in{1,...,M} for sa and then computing \prod{j=1,...,k-1}{1-h_j(sW))*h_k(sW)}
  #    which is equivalent to the joint conditional probability that sa belongs to the interval (i_k,i_{k+1}) and does not belong
  #    to any of the intervals 1 to k-1, conditional on sW. 
  #    The evaluation above utilizes a discretization of the fact that any continous density f of random variable X can be written as f_X(x)=S_X(x)*h_X(x), 
  #    for a continuous density f of X where S_X(x):=P(X>x) is the survival function for X, h_X=P(X>x|X>=x) is the hazard function for X; as well as the fact that
  #    the discretized survival function S_X(x) can be written as a of the hazards for s<x: S_X(x)=\prod{s<x}h_X(x).
## **********************************************************************
## Two methods for discretizing (creating bin intervals) for a continuous summary measure sA
## **********************************************************************
# There are 2 alternative methods to defining the bin cutoffs S=(i_1,...,i_M,i_{M+1}) for a continuous summary measure sA.
# The choice of which method is used along with other discretization parameters (e.g., total number of bins) is controlled via the tmlenet_options() function. 
# See ?tmlenet_options for additional details.
  # *********************
  # Approach 1 (default, equal length intervals):
  # *********************
  # The bins are defined by splitting the range of observed sA (sa_1,...,sa_n) into equal length intervals. 
  # This is the dafault discretization method, set by passing an argument \code{binByMass=FALSE} to the \code{tmlenet_options} function.
  # The intervals will be defined by splitting the range of (sa_1,...,sa_N) into \code{nbins} number of equal length intervals, 
  # where \code{nbins} is another argument of \code{tmlenet_options()} function.
  # When \code{nbins=NA} (the default setting) the actual value of \code{nbins} is computed at run time by taking the integer value (floor) of \code{n/maxNperBin},
  # for \code{n} - the total observed sample size and \code{maxNperBin=1000} - another argument of \code{tmlenet_options()} with the default value 1,000.
  # *********************
  # Approach 2 (data-adaptive unequal length intervals):
  # *********************
  # The intervals are defined by splitting the range of \code{sA} into non-equal length data-adaptive intervals that ensures that each interval contains around 
  # \code{maxNperBin} observations from (sa_j:j=1,...,N).
  # This interval definition approach is activated by passing an argument \code{binByMass = TRUE} to \code{tmlenet_options()}.
  # The method ensures that an approximately equal number of observations will belong to each interval, where that number of observations for each interval
  # is controlled by setting \code{maxNperBin}. The default setting is \code{maxNperBin=1000} observations per interval.
  # *********************
  # Approach 3 (data-adaptive dhist approach that is a mix of Approaches 1 & 2)
  # *********************
## ---------------------------------------------------------------------

# **********
# TO DO (ContinOutModel)
# x) See how to generalize to pooled fits, k-specific fits, etc (use subset definitions + ?)
# **********
# TO DO (ContinOutModel, binirize)
# * (***BUG***) Currently make.bins_mtx_1 fails on binary A with automatic bin detection (all values are placed in last (2nd) bin)
# * Implement regression for categorical outvar: the process is identical to contin, except that normalize, define.intervals() & discretize() is skipped
# * Need to test that make.bins_mtx_1 will do the right thing when x_cat is (0, ..., ncats) instead of (1, ..., ncats) (IT SHOULD WORK)

# Generic S3 constructor for the summary model classes:
NewSummaryModel <- function(reg, O.datnetA, ...) { UseMethod("NewSummaryModel") }
# Summary model constructor for binary outcome sA[j]:
NewSummaryModel.binary <- function(reg, ...) {
  if (gvars$verbose) print("BinOutModel constructor called...")
  BinOutModel$new(reg = reg, ...)
}
# Summary model constructor for continuous outcome sA[j]:
NewSummaryModel.contin <- function(reg, O.datnetA, ...) { 
  if (gvars$verbose) print("ContinOutModel constructor called...")
  ContinOutModel$new(reg = reg, O.datnetA = O.datnetA, ...)
}

RegressionClass <- R6Class("RegressionClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    outvar.class = character(),    # vector for classes of the outcome vars: bin / cont / cat
    outvar = character(),          # vector of regression outcome variable names
    predvars = character(),        # either a pool of all predictors (sW) or regression-specific predictor names
    reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard)
    subset = NULL,                 # subset expression (later evaluated to logical vector in the envir of the data)
    ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    nbins = NULL,                  # actual nbins used, for cont. outvar, defined in ContinOutModel$new()
    bin_nms = NULL,                # column names for bin indicators
    useglm = logical(),            # TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
    parfit = logical(),            # TRUE for fitting binary regressions in parallel
    bin_bymass = logical(),        # for cont outvar, create bin cutoffs based on equal mass distribution?
    bin_bydhist = logical(),       # if TRUE, use dhist approach for bin definitions
    max_nperbin = integer(),       # maximum n observations allowed per binary bin
    pool_cont = logical(),         # Pool binned cont outvar obs into long format (adding bin_ID as a covaraite)
    outvars_to_pool = character(), # Names of the binned continuous sVars, should match bin_nms

    # NAMED VECTOR THAT CONTAINS (bw_j : j=1,...,M). CAN BE QUERIED BY BinOutModel$predictAeqa() as: intrvls.width[outvar]
    # INSIDE BinOutModel$predictAeqa() EACH probA1 is adjusted as probA1/self$reg$intrvls.width[self$getoutvarnm]
    # WHERE outvar is a discretized bin name for continuous sA
    # BinOutModel$predictAeqa() THIS CLASS
    # When sA is not continuous, intrvls.width IS SET TO 1.
    # When sA is continuous, intrvls.width is SET TO self$intrvls.width INSIDE ContinOutModel$new() with names(intrvls.width) <- reg$bin_nms
    intrvls.width = 1L,
    intrvls = numeric(),
    # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
    # form = NULL,                 # (NOT IMPLEMENTED) reg formula, if provided run using the usual glm / speedglm functions
    initialize = function(outvar.class = gvars$sVartypes$bin,
                          outvar, predvars, subset,
                          # Needed to add ReplMisVal0 = TRUE for case sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise:
                          ReplMisVal0 = TRUE,  
                          useglm = getopt("useglm"),
                          parfit = getopt("parfit"),
                          nbins = getopt("nbins"),
                          bin_bymass = getopt("binByMass"),
                          bin_bydhist = getopt("binBydhist"),
                          max_nperbin = getopt("maxNperBin"),
                          pool_cont = getopt("poolContinVar")
                          ) {

      self$outvar.class <- outvar.class
      self$outvar <- outvar
      self$predvars <- predvars
      self$ReplMisVal0 <- ReplMisVal0

      self$useglm <- useglm
      self$parfit <- parfit

      self$nbins <- nbins
      self$bin_bymass <- bin_bymass
      self$bin_bydhist <- bin_bydhist
      self$max_nperbin <- max_nperbin
      self$pool_cont <- pool_cont

      n_regs <- length(outvar)
      if (!missing(subset)) {
        self$subset <- subset
        if (length(subset) < n_regs) {
          self$subset <- rep_len(subset, n_regs)
        } else if (length(subset) > n_regs) {
          # ... TO FINISH ...
          if (!is.logical(subset)) stop("not implemented")
          if (gvars$verbose) message("logical subset index length: " %+% length(subset))
          # increase n_regs to all combinations of (n_regs x subset)
        }
      } else {
        self$subset <- rep_len(list(TRUE), n_regs)
      }
    },

    setToKthRegresssion = function(k_i, reg, regs_list) {
      if (!missing(k_i)) {
        if (missing(reg)) stop("reg must be also specified when k_i is specified")
        assert_that(k_i <= length(reg$outvar))  
        self$outvar.class <- reg$outvar.class[[k_i]] # Class of the outcome var: binary, categorical, continuous:
        self$outvar <- reg$outvar[[k_i]] # An outcome variable that is being modeled:
        n_regs <- length(reg$outvar)

        if (self$reg_hazard) {
          covars_nms <- reg$predvars # Regression covars (predictors):  
        } else {
          covars_nms <- c(reg$outvar[-c(k_i:n_regs)], reg$predvars) # Regression covars (predictors):  
        }
        self$predvars <- covars_nms

        if (is.list(self$subset)) {
          self$subset <- reg$subset[[k_i]]
        } else {
          self$subset <- reg$subset
        }
        self$S3class <- self$outvar.class # Set class on self for S3 dispatch...
      } else {
        self$outvar.class <- regs_list$outvar.class # Vector of class(es) of outcome var(s): binary, categorical, continuous
        self$outvar <- regs_list$outvar # An outcome variable that is being modeled:
        self$predvars <- regs_list$predvars
        self$subset <- regs_list$subset
      }
    },

    resetS3class = function() class(self) <- c("RegressionClass", "R6")

  ),

  active = list(
    # For S3 dispatch on NewSummaryModel():
    S3class = function(newclass) {
      if (!missing(newclass)) {
        if (length(newclass) > 1) stop("cannot set S3 class on RegressionClass with more than one outvar variable")
        if (length(class(self)) > 2) stop("S3 dispatch class on RegressionClass has already been set")
        class(self) <- c(class(self), newclass)
      } else {
        return(class(self))
      }
    },

    get.reg = function() {
      list(outvar.class = self$outvar.class,
          outvar = self$outvar,
          predvars = self$predvars,
          subset = self$subset)
    }
  )
)

## ---------------------------------------------------------------------
# Class for defining, fitting and predicting the probability model P(sA = sa | sW = sw) under g.star or under g.0 for summary measures (sW,sA)
# Accepts (1) data (data.frame) for (sA,sW), (2) newdata (data.frame) for prediction, (3) obsdat.sA (matrix) for sa values;
# Defines and manages the factorization of the joint P(sA = sa | ... ) into reg models sA[j] ~ \bar{sA[j-1]} + sW;
# Figures out reg mdel factorization based on name ordering in (sA_nms, sW_nms);
# Evaluates subset_exprs in the envirs of data and newdata data.frames
# Calls BinOutModel$new, assumes each sA[j] is binary in reg (sA[j] ~ \bar{sA[j-1]} + sW);
## ---------------------------------------------------------------------
#' @title R6 class for fitting and predicting model P(sA|sW) under g.star or g.0
#' @docType class
#' @format An R6 class object.
#' @name SummariesModel
#' @details Class for defining, fitting and predicting the probability model P(sA = sa | sW = sw) under g_star or under g_0 for summary measures (sW,sA).
#' \itemize{
#' \item{n_regs} ...
#' \item{sA_nms} ...
#' \item{sW_nms} ...
#' }
#' Additional details about implementation of the class...
# @export
SummariesModel <- R6Class(classname = "SummariesModel",
	portable = TRUE,
	class = TRUE,
	public = list(
		n_regs = integer(),        # total no. of reg. models (logistic regressions)
    parfit_allowed = FALSE,    # allow parallel fit of multivar outvar when 1) reg$parfit = TRUE & 2) all.outvar.bin = TRUE
    initialize = function(reg, ...) {
			self$n_regs <- length(reg$outvar) # Number of sep. logistic regressions to run
      all.outvar.bin <-  all(reg$outvar.class %in% gvars$sVartypes$bin)

      if (reg$parfit & all.outvar.bin & (self$n_regs > 1)) self$parfit_allowed <- TRUE

      # if (gvars$verbose) {
        print("#----------------------------------------------------------------------------------");
        print("New SummariesModel object:");
        print("No. of regressions: " %+% self$n_regs)
        # print("sA_classes (reg$outvar.class): " %+% paste(reg$outvar.class, collapse = ", "))
        print("All outvar binary? " %+% all.outvar.bin)
        print("Doing parallel fit? " %+% self$parfit_allowed)
        print("sW_nms (reg$predvars): " %+% paste(reg$predvars, collapse = ", "))
        print("#----------------------------------------------------------------------------------");        
      # }

      # factorize the joint into univariate regressions, by dimensionality of the outcome variable (sA_nms):
			for (k_i in 1:self$n_regs) {
        reg_i <- reg$clone()
        reg_i$setToKthRegresssion(k_i, reg)
        # Calling the constructor for the summary model P(sA[j]|\bar{sA}[j-1], sW}), dispatching on reg_i class
        PsAsW.model <- NewSummaryModel(reg = reg_i, ...)
				private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
				names(private$PsAsW.models)[k_i] <- "P(sA|sW)."%+%k_i
			}
			invisible(self)
		},

		length = function(){ base::length(private$PsAsW.models) },
		getPsAsW.models = function() { private$PsAsW.models },  # get all summary model objects (one model object per outcome var sA[j])
		getcumprodAeqa = function() { private$cumprodAeqa },  # get joint prob as a vector of the cumulative prod over j for P(sA[j]=a[j]|sW)

    # **********************************************************************
    # TO DO: Add $copy.fit(SummariesModel) method that will propagate copy all the model fits down the line
    copy.fit = function(SummariesModel) {},
    # **********************************************************************

    fit = function(data) {
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$fit(data = data)
        }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each fitRes[[k_i]] will contain a copy of every single R6 object that was passed by reference -> 
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        fitRes <- foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$fit(data = data)
        }
        # copy the fits one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.fit(fitRes[[k_i]])
        }
      }
		  invisible(self)
		},

    # TO DO rename to:
    # predictP_1 = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
		predict = function(newdata) {
		  if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
		    # return(invisible(self))
        stop("must provide newdata")
		  }
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
		    for (k_i in seq_along(private$PsAsW.models)) { 
		      private$PsAsW.models[[k_i]]$predict(newdata = newdata)
		    }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each predRes[[k_i]] will contain a copy of every single R6 object that was passed by reference -> 
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        predRes <- foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        }
        # copy the predictions one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.predict(predRes[[k_i]])
        }
        # private$PsAsW.models[] <- predRes
      }
		  invisible(self)
		},


		# WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
		# Uses daughter objects (stored from prev call to fit()) to get predictions for P(sA=obsdat.sA|sW=sw)
		# Invisibly returns the joint probability P(sA=sa|sW=sw), also saves it as a private field "cumprodAeqa"
		predictAeqa = function(obs.DatNet.sWsA, ...) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
			assert_that(!missing(obs.DatNet.sWsA))
      assert_that(is.DatNet.sWsA(obs.DatNet.sWsA))
			n <- obs.DatNet.sWsA$nobs

      if (!self$parfit_allowed) {
			 cumprodAeqa <- rep.int(1, n)
			 for (k_i in seq_along(private$PsAsW.models)) { # loop over all regressions in PsAsW.models
				    cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA, ...)
			 }
      } else if (self$parfit_allowed) {
        mcoptions <- list(preschedule = TRUE)
        probAeqa_list <- foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA, ...)
        }
        # print("probAeqa_list: "); print(class(probAeqa_list)); print(length(probAeqa_list)); print(head(probAeqa_list[[1]]))
        # probAeqa <- data.table::rbindlist(probAeqa)
        cbind_t <- system.time(
          probAeqa_mat <- do.call('cbind', probAeqa_list)
          )
        print("cbind_t: "); print(cbind_t)
        # print("probAeqa_mat: "); print(dim(probAeqa_mat)); print(head(probAeqa_mat)); print(class(probAeqa_mat))
        rowProds_t <- system.time(
          cumprodAeqa <- matrixStats::rowProds(probAeqa_mat)
          )
        print("rowProds_t: "); print(rowProds_t)
        # print("cumprodAeqa"); print(head(cumprodAeqa)); print(length(cumprodAeqa))
      }
      private$cumprodAeqa <- cumprodAeqa
			return(cumprodAeqa)
		}
	),

	active = list(
    # recursively call all saved daughter model fits and wipe out any traces of saved data
    wipe.alldat = function() {  
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$wipe.alldat
      }
      return(self)
    },
		actplhold = function() {} # placeholder, not used
	),

	private = list(
		PsAsW.models = list(),
		fitted.pbins = list(),
		cumprodAeqa = NULL
	)
)

## ---------------------------------------------------------------------
# R6 clas inherits from SummariesModel, fitting for continuous outcome sA[j] in sA
# Called from SummariesModel for contin sA[j]. Gets passed new subset definitions, e.g., (!mis(Bin_sA[j]_i))
# Defines the fitting algorithm for probability sA[j] ~ \bar{sA[j-1]} + sW
# Reconstructs the likelihood P(sA[j]=sa[j]|sW)
# Continuous sA[j] is discretized using data-adaptive bin defintions
# then estimates regressions for the hazard Bin_sA[j][i] ~ \bar{sA[j-1]} + sW.
# I.e., we estimate the probability that continuous sA[j] falls into bin Bin_sA[j]_i,
# given that sA[j] is not in bins Bin_sA[1]_1, ..., Bin_sA[1]_{i-1}.
# The dataset of discretized summary measures (BinsA[j]_1,...,BinsA[j]_M) is created 
# while discretizing sA[j] into M bins.
## ---------------------------------------------------------------------
ContinOutModel <- R6Class(classname = "ContinOutModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    # O.datnetA = NULL,
    reg = NULL,
    outvar = character(),     # the name of the continous outcome var (sA[j])
    intrvls = NULL,
    intrvls.width = NULL,
    bin_weights = NULL,

    # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
    initialize = function(reg, O.datnetA, datnet.gstar, ...) {
      # self$O.datnetA <- O.datnetA
      self$reg <- reg
      self$outvar <- reg$outvar
      self$intrvls <- O.datnetA$detect.sVar.intrvls(reg$outvar,
                                                    nbins = self$reg$nbins,
                                                    bin_bymass = self$reg$bin_bymass,
                                                    bin_bydhist = self$reg$bin_bydhist,
                                                    max_nperbin = self$reg$max_nperbin)
      if (!missing(datnet.gstar)) {
        message("Found datnet.gstar! sVar intervals under gstar will be based on a union of intervals under g0 and gstar")
        gstar.intrvls <- datnet.gstar$detect.sVar.intrvls(reg$outvar,
                                                    nbins = self$reg$nbins,
                                                    bin_bymass = self$reg$bin_bymass,
                                                    bin_bydhist = self$reg$bin_bydhist,
                                                    max_nperbin = self$reg$max_nperbin)
        self$intrvls <- unique(sort(union(self$intrvls, gstar.intrvls)))
      }

      # Define the number of bins (no. of binary regressions to run), new outvar var names (bin names) # all predvars remain unchanged
      self$reg$intrvls <- self$intrvls
      self$reg$nbins <- length(self$intrvls) - 1
      self$reg$bin_nms <- O.datnetA$bin.nms.sVar(reg$outvar, self$reg$nbins)

      # Save bin widths in reg class (naming the vector entries by bin names):
      self$intrvls.width <- diff(self$intrvls)
      self$intrvls.width[self$intrvls.width <= gvars$tolerr] <- 1
      self$reg$intrvls.width <- self$intrvls.width
      names(self$reg$intrvls.width) <- names(self$intrvls.width) <- self$reg$bin_nms

      # INSTEAD OF DEFINING NEW RegressionClass now cloning parent reg object and then ADDING new SETTINGS:
      bin_regs <- self$reg$clone()
      bin_regs$reg_hazard <- TRUE # to not add previous degenerate bin indicators as predictor covariates for each bin regression

      # if (gvars$verbose)  {
        print("ContinOutModel sA: "%+%self$outvar)
        # print("ContinOutModel bin intervals:"); print(self$intrvls)
        print("ContinOutModel reg$nbins: " %+% self$reg$nbins)
      # }

      # Define subset evaluation for new bins:
      if (!self$reg$pool_cont) {
        add.oldsubset <- TRUE
        new.subsets <- lapply(self$reg$bin_nms,
                                  function(var) { 
                                    res <- var
                                    if (add.oldsubset) res <- c(res, self$reg$subset)
                                    res
                                  })

        new.sAclass <- as.list(rep_len(gvars$sVartypes$bin, self$reg$nbins))
        names(new.sAclass) <- self$reg$bin_nms
        bin_regs$setToKthRegresssion(regs_list = list(outvar.class = new.sAclass, 
                                                      outvar = self$reg$bin_nms, 
                                                      predvars = self$reg$predvars, 
                                                      subset = new.subsets))

      } else {
        bin_regs$outvar.class <- gvars$sVartypes$bin
        bin_regs$outvar <- self$outvar
        bin_regs$outvars_to_pool <- self$reg$bin_nms

        if (gvars$verbose)  {
          print("pooled bin_regs$outvar: "); print(bin_regs$outvar)
          print("bin_regs$outvars_to_pool: "); print(bin_regs$outvars_to_pool)
          print("bin_regs$subset: "); print(bin_regs$subset)
        }
      }

      bin_regs$resetS3class()
      super$initialize(reg = bin_regs, ...)
    },

    # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      if (gvars$verbose) message("fit in continuous summary... for: " %+% self$outvar)
      # input data is of class DatNet.sWsA
      # data$datnetW$dat.sVar -> interface to get sW (currently matrix)
      # data$datnetA$dat.sVar -> interface to get sA (currently matrix)
      # data$dat.sVar -> returns a combined mat of sWsA (no subsetting, no covar sel)
      # data$get.dat.sWsA(rowsubset, covars) -> returns a processed mat of sWsA

      if (gvars$verbose) print("current active bin sVar: " %+% data$active.bin.sVar)
      # Binirizes & saves binned matrix inside DatNet.sWsA
      data$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
 
      if (gvars$verbose) {
        print("active bin sVar after calling binirize.sVar: " %+% data$active.bin.sVar)
        print("data$dat.sVar for: " %+% self$outvar); print(head(data$dat.sVar, 5))
        print("binned dataset for: " %+% self$outvar); print(head(cbind(data$ord.sVar, data$dat.bin.sVar), 5))
      }
      print("freq count for transformed ord.sVar: "); print(table(data$ord.sVar))

      super$fit(data) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")

      data$emptydat.bin.sVar # wiping out binirized mat in data DatNet.sWsA object...
      self$wipe.alldat # wiping out all data traces in ContinOutModel...
      invisible(self)
    },

    # TO DO: rename to:
    # predictP_1 = function(newdata, response = FALSE) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata, response = FALSE) {
      if (missing(newdata)) {
        stop("must provide newdata")
        # return(invisible(self))
      }
      if (gvars$verbose) message("predict in continuous summary...")
      # mat_bin doesn't need to be saved (even though its invisibly returned); mat_bin is automatically saved in datnet.sW.sA - a potentially dangerous side-effect!!!
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      # if (gvars$verbose) {
      #   print("mat_bin"); print(head(mat_bin))
      # }
      super$predict(newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DatNet.sWsA object...
      invisible(self)      
    },

    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(obs.DatNet.sWsA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      if (gvars$verbose) message("predictAeqa in continuous summary...")
      assert_that(is.DatNet.sWsA(obs.DatNet.sWsA))
      # print("current active bin sVar: " %+% obs.DatNet.sWsA$active.bin.sVar)
      obs.DatNet.sWsA$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      bws <- obs.DatNet.sWsA$get.sVar.bw(name.sVar = self$outvar, intervals = self$intrvls)
      self$bin_weights <- (1 / bws) # weight based on 1 / (sVar bin widths)
      # print("mat_bin: "); print(dim(mat_bin))
      # print(head(cbind(obs.DatNet.sWsA$get.outvar(var = self$outvar), bws = bws, bw_diff = bw_diff)))
      # print("bws"); print(length(bws))

      # Option 1: ADJUST FINAL PROB by bw.j TO OBTAIN density at a point f(sa|sw) = P(sA=sa|sW=sw):
      cumprodAeqa <- super$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA) * self$bin_weights
      # Alternative 2: ALso integrate the difference of sA value and its left most bin cutoff: x - b_{j-1} and pass it
      # This is done so that we can integrate the constant hazard all the way to the value of x:
        # * (1 - bw.j.sA_diff*(1/self$bin_weights)*probA1) (discrete)
        # * exp(-bw.j.sA_diff*(1/self$bin_weights)*probA1) (continuous)
      # bw.j.sA_diff <- obs.DatNet.sWsA$get.sVar.bwdiff(name.sVar = self$outvar, intervals = self$intrvls)
      # cumprodAeqa <- super$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA, bw.j.sA_diff = bw.j.sA_diff) * self$bin_weights

      obs.DatNet.sWsA$emptydat.bin.sVar # wiping out binirized mat in obs.DatNet.sWsA object...
      self$bin_weights <- NULL # wiping out self$bin_weights...
      self$wipe.alldat # wiping out all data traces in ContinOutModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)

