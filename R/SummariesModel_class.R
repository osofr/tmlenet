#' @importFrom assertthat assert_that

# **********
# TO DO (ContinOutModel)
# x) Remove subset_expr definition from here, where is a more approapriate location?
# x) See how to generalize to pooled fits, k-specific fits, etc (use subset definitions + ?)
# x) Need a function to convert data_mtx to long format for SummaryM.pool
# x) For SummaryM.pool how to define predict function?

# **********
# TO DO (ContinOutModel, binirize)
# * (***BUG***) Currently make.bins_mtx_1 fails on binary A with automatic bin detection (all values are placed in last (2nd) bin)
# * Implement regression for categorical outvar: the process is identical to contin, except that normalize, define.intervals() & discretize() is skipped
# * Need to test that make.bins_mtx_1 will do the right thing when x_cat is (0, ..., ncats) instead of (1, ..., ncats) (IT SHOULD WORK)
# * Create a class RegressionContClass that inherits from RegressionClass for continous outcome sVar?

RegressionClass <- R6Class("RegressionClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    useglm = FALSE,             # TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
    outvar.class = character(), # vector for classes of the outcome vars: bin / cont / cat
    outvar = character(),       # vector of regression outcome variable names
    predvars = character(),     # either a pool of all predictors (sW) or regression-specific predictor names
    subset = NULL,              # subset expression (later evaluated to logical vector in the envir of the data)
    ReplMisVal0 = TRUE,         # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    nbins = NULL,               # for cont. outvar, copied by ContinOutModel$new from O.datnetA$nbins.sVar(reg$outvar)
    bin_nms = NULL,             # column names for bin indicators
    bin_bymass = TRUE,          # for cont outvar create bins based on equal interval cutoffs rather than equal mass cutoffs
    max_nperbin = gvars$max_nperbin,
    bin_width = 0L,
    pool_cont = FALSE,          # (NOT IMPLEMENTED) pool binned cont outvar obs into long format (adding bin_n as a covaraite)
    family = NULL,              # (NOT IMPLEMENTED) to run w/ other than "binomial" family
    # form = NULL,              # (NOT IMPLEMENTED) reg formula, if provided run using the usual glm / speedglm functions
    initialize = function(outvar.class = gvars$sVartypes$bin,
                          outvar, predvars, subset, useglm = FALSE,
                          ReplMisVal0 = TRUE,  # Added ReplMisVal0 = TRUE for cases such as sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise:
                          bin_bymass = TRUE,
                          max_nperbin = gvars$max_nperbin) {

      self$outvar.class <- outvar.class
      self$outvar <- outvar
      self$predvars <- predvars
      self$useglm <- useglm
      self$ReplMisVal0 <- ReplMisVal0
      self$bin_bymass <- bin_bymass
      self$max_nperbin <- max_nperbin

      n_regs <- length(outvar)
      if (!missing(subset)) {
        self$subset <- subset
        if (length(subset) < n_regs) {
          self$subset <- rep_len(subset, n_regs)
        } else if (length(subset) > n_regs) {
          # ... TO FINISH ...
          if (!is.logical(subset)) stop("not implemented")
          message("logical subset index length: " %+% length(subset))
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
        covars_nms <- c(reg$outvar[-c(k_i:n_regs)], reg$predvars) # Regression covars (predictors):
        self$predvars <- covars_nms
        self$subset <- reg$subset[[k_i]]
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

# Template for RegressionContClass (continous outvar):
# NOT IMPLEMENTED YET
# RegressionContClass <- R6Class("RegressionContClass",
#   inherit = RegressionClass,
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     nbins = NULL,               # for cont. outvar, copied by ContinOutModel$new from O.datnetA$nbins.sVar(reg$outvar)
#     bin_nms = NULL,             # column names for bin indicators
#     bin_bymass = FALSE,         # for cont outvar create bins based on equal interval cutoffs rather than equal mass cutoffs
#     pool_cont = FALSE,          # (NOT IMPLEMENTED) pool binned cont outvar obs into long format (adding bin_n as a covaraite)
#     initialize = function(bin_bymass = FALSE, pool_cont = TRUE, ...) {
#       self$bin_bymass <- bin_bymass
#       self$pool_cont <- pool_cont
#       super$initialize(...) # call the parent class contructor
#     }
#   ),
#   active = list(
#     plchlder = function() { } # placeholder, not used
#   )
# )

# Generic S3 constructor for the summary model classes:
NewSummaryModel <- function(reg, O.datnetA, ...) { UseMethod("NewSummaryModel") }
# Summary model constructor for binary outcome sA[j]:
NewSummaryModel.binary <- function(reg, ...) {
  # print("BinOutModel constructor called...")
  # BinOutModel$new(glm = TRUE, reg = reg, ...) # fit a model with new object BinOutModel using glm.fit
  BinOutModel$new(glm = FALSE, reg = reg, ...) # fit a model with new object BinOutModel using speedglm.wfit
  # Alternative name: BinOutModel
}
# Summary model constructor for categorical outcome sA[j]:
# NewSummaryModel.cat = function(reg, ...) { 
#   print("CatSummaryModel constructor called...")
#   CatSummaryModel$new(...)
# }
# Summary model constructor for continuous outcome sA[j]:
NewSummaryModel.contin <- function(reg, O.datnetA, ...) { 
	print("ContinOutModel constructor called...")
	ContinOutModel$new(reg = reg, O.datnetA = O.datnetA, ...)
	# Alternative name: ContinOutModel
}

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
		n_regs = integer(),     # total no. of reg. models (logistic regressions)
    initialize = function(reg, ...) {
			self$n_regs <- length(reg$outvar) # Number of sep. logistic regressions to run
			print("#----------------------------------------------------------------------------------");
			print("New SummariesModel object:");
			print("No. of regressions: " %+% self$n_regs)
			print("sA_nms (reg$outvar): " %+% paste(reg$outvar, collapse = ", "))
			print("sW_nms (reg$predvars): " %+% paste(reg$predvars, collapse = ", "))
			print("#----------------------------------------------------------------------------------");

      # factorize by dimensionality of the outcome variable (sA_nms):
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

    fit = function(data) {
		  for (k_i in seq_along(private$PsAsW.models)) { # loop over all regressions in PsAsW.models
		    private$PsAsW.models[[k_i]]$fit(data = data) # below is replaced with this
		  }
		  invisible(self)
		},

		predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
		  if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
		    return(invisible(self))
		  }
		  for (k_i in seq_along(private$PsAsW.models)) { # loop over all regressions in PsAsW.models
		    private$PsAsW.models[[k_i]]$predict(newdata = newdata)
		  }
		  invisible(self)
		},

		# WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
		# Uses daughter objects (stored from prev call to fit()) to get predictions for P(sA=obsdat.sA|sW=sw)
		# Invisibly returns the joint probability P(sA=sa|sW=sw), also saves it as a private field "cumprodAeqa"
		predictAeqa = function(obs.DatNet.sWsA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
			assert_that(is.DatNet.sWsA(obs.DatNet.sWsA))
			n <- obs.DatNet.sWsA$nobs
			cumprodAeqa <- rep_len(1, n)
			for (k_i in seq_along(private$PsAsW.models)) { # loop over all regressions in PsAsW.models
				cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA)
			}
			private$cumprodAeqa <- cumprodAeqa
			invisible(cumprodAeqa)
		}
	),

	active = list(
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
# Defines the fitting algorithm for sA[j] ~ \bar{sA[j-1]} + sW
# Reconstructs the likelihood P(sA[j]=sa[j]|sW)
# If needed, discretizes the continuous sA[j] and correctly defines regressions BinsA[j][i] ~ \bar{BinsA[j][i-1]} + \bar{sA[j-1]} + sW
# Creates the appropriate dataset of discretized summary measures (BinsA[j][1],...,BinsA[j][M])
# Called from SummariesModel for contin sA[j]. Gets passed new subset definitions, e.g., (!mis(BinsA))
## ---------------------------------------------------------------------
ContinOutModel <- R6Class(classname = "ContinOutModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    O.datnetA = NULL,
    reg = NULL,
    outvar = character(),     # the name of the continous outcome var (sA[j])
    intrvls = NULL,
    intrvls.width = NULL,
    bin_weights = NULL,

    # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
    initialize = function(reg, O.datnetA, ...) {
      # Define the number of bins (no. of regs to run), new outvar var names, predvars remain unchanged:
      self$O.datnetA <- O.datnetA
      self$reg <- reg
      self$outvar <- reg$outvar
      self$reg$nbins <- O.datnetA$nbins.sVar(reg$outvar)

      if (self$reg$nbins < 2) {
        # stop("bin intervals were not defined for continuous sA: " %+% reg$outvar)
        message("defining bin intervals that were not defined for continuous sA: " %+% reg$outvar)
        self$intrvls <- O.datnetA$detect.sVar.intrvls(reg$outvar, 
                                                      bin_bymass = self$reg$bin_bymass, 
                                                      max_nperbin = self$reg$max_nperbin)

        self$intrvls.width <- diff(self$intrvls)
        print("defined bin intervals:"); print(self$intrvls)
        O.datnetA$set.sVar.intrvls(reg$outvar, self$intrvls)
        self$reg$nbins <- O.datnetA$nbins.sVar(reg$outvar)
      }

      self$reg$bin_nms <- O.datnetA$bin.nms.sVar(reg$outvar) # new.sA_nms <- O.datnetA$bin.nms.sVar(reg$outvar)
      new.sAclass <- as.list(rep_len(gvars$sVartypes$bin, self$reg$nbins))
      names(new.sAclass) <- self$reg$bin_nms


      # #todo 27 (ContinOutModel) +0: Put subset eval in a separate function (with var as arg + additional args)
      add.oldsubset <- TRUE
      new.subsets_chr <- lapply(self$reg$bin_nms,
                                function(var) {
                                  newsub <- "!misfun("%+%var%+%")"
                                    if (add.oldsubset) {
                                      newsub <- newsub %+% " & (" %+% deparse(reg$subset) %+% ")"
                                    }
                                    newsub
                                })
      new.subsets <- lapply(new.subsets_chr,
                                function(subset_chr) {
                                  subset_expr <- try(parse(text=subset_chr)[[1]])
                                  if(inherits(subset_expr, "try-error")) stop("can't parse the subset formula", call.=FALSE)
                                  subset_expr
                                })

      # INSTEAD OF DEFINING NEW RegressionClass now cloning reg and then COPY new SETTINGS:
      bin_regs <- reg$clone()
      bin_regs$setToKthRegresssion(regs_list = list(outvar.class = new.sAclass, 
                                                    outvar = self$reg$bin_nms, 
                                                    predvars = reg$predvars, 
                                                    subset = new.subsets))
      bin_regs$resetS3class()

      print("contin sA: "%+%self$outvar);
      print("contin sA reg$nbins: " %+% self$reg$nbins);
      print("contin binned sA names: " %+% paste(self$reg$bin_nms, collapse = ",")); 
      super$initialize(reg = bin_regs, ...)
    },

    # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      message("fit in continuous summary... for: " %+% self$outvar)
      # input data is of class DatNet.sWsA
      # data$datnetW$dat.sVar -> interface to get sW (currently matrix)
      # data$datnetA$dat.sVar -> interface to get sA (currently matrix)
      # data$dat.sVar -> returns a combined df / mat? of sWsA (no subsetting, no covar sel)
      # data$get.dat.sWsA(rowsubset, covars) -> returns a processed df / mat? of sWsA

      # NOTE, THIS LINE IS VERY IMPORTANT
      # IF bins were defined inside ContinOutModel$new(), they will be saved in O.datnetA and DatNet.sWsA will not see them
      # unless the intervals are copied from O.datnetA
      data$copy.cbin.intrvls()
      print("intervals after data$copy.cbin.intrvls() in ContinOutModel$fit: "); print(data$get.sVar.intrvls(self$outvar))
      print("current active bin sVar: " %+% data$active.bin.sVar)
      # Saves binned matrix to DatNet.sWsA:
      data$binirize.sVar(self$outvar) # Note, don't need to save mat_bin here, its already saved inside DatNet.sWsA: # mat_bin <- data$binirize.sVar(self$outvar)
      print("active bin sVar after calling binirize.sVar: " %+% data$active.bin.sVar)
      ord.sVar <- data$discretize.sVar(self$outvar)
      print("data$dat.sVar for: " %+% self$outvar); print(head(data$dat.sVar, 5))
      print("binned dataset for: " %+% self$outvar); print(head(cbind(ord.sVar, data$binirize.sVar(self$outvar)), 5))
      # print("freq count for original variable: "); print(table(data$dat.sVar[, self$outvar]))
      print("freq count for transformed ord.sVar: "); print(table(ord.sVar))

      super$fit(data) # call the parent class fit method
      message("fit for " %+% self$outvar %+% " var succeeded...")
      invisible(self)
    },

    predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
      if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
        return(invisible(self))
      }
      print("predict in continuous summary...")

      newdata$copy.cbin.intrvls() #07/18/15 added copy of interval definitions:

      # NEW VERSION. mat_bin doesn't need to be saved (even though its invisibly returned)
      # mat_bin is automatically saved in datnet.sW.sA - potentially dangerous due to side-effects!!!
      mat_bin <- newdata$binirize.sVar(self$outvar)
      print("mat_bin"); print(head(mat_bin))
      super$predict(newdata)
      invisible(self)
    },

    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(obs.DatNet.sWsA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      print("predictAeqa in continuous summary...")
      # added 07/18/15:
      obs.DatNet.sWsA$copy.cbin.intrvls()

      assert_that(is.DatNet.sWsA(obs.DatNet.sWsA))
      assert_that(is.vector(obs.DatNet.sWsA$get.outvar(var = self$outvar)))
      
      obs.DatNet.sWsA$binirize.sVar(self$outvar)
      # mat_bin <- obs.DatNet.sWsA$binirize.sVar(self$outvar)
      bws <- obs.DatNet.sWsA$get.sVar.bw(self$outvar)
      self$bin_weights <- (1 / bws) # weight based on 1 / (sVar bin widths)
      cumprodAeqa <- super$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA) * self$bin_weights
      invisible(cumprodAeqa)
      # invisible(super$predictAeqa(self$binirize(obsdat.sA)))  # one line alternative to the above?
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)

## ---------------------------------------------------------------------
# NOT IMPLEMENTED
# Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] 
# add a covariate BinInd (1,..,M), convert to long format (turn BinsA[1], ...., BinsA[M] into rows) 
# then call $super$fit on transformed data?
## ---------------------------------------------------------------------
ContinOutModelPool <- R6Class(classname = "ContinOutModelPool",
  inherit = ContinOutModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
    initialize = function(...) {
    	super$initialize(...) # call the parent class contructor
    }
  ),
  active = list(
    actplhold = function() {} # placeholder
  ),
  private = list(
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)

# THIS IS A BAD APPROACH (DEFINING SEPARATE CLASS FOR CAT). NEED TO HANDLE cat sA[j] INSIDE ContinOutModel
# CatSummaryModel <- R6Class(classname = "CatSummaryModel",
#   inherit = SummariesModel,
#   portable = TRUE,
#   class = TRUE
# )
