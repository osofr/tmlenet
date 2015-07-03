# **********
# TO DO (RegressionClass)
# x) It would be much nicer for SummariesModel$new() to accept only (reg, ...) like the rest of Model constructors.
# x) RegressionClass could be also defined to have vectors for outvar.class, outvar, subset
# x) Will then be able to use glm flag in RegressionClass as the only control tool for speedglm vs. glm (currently this flag is unused)

# **********
# TO DO (ContinSummaryModel)
# x) Remove subset_expr definition from here, where is a more approapriate location?
# x) See how to generalize to pooled fits, k-specific fits, etc (use subset definitions + ?)
# x) Need a function to convert data_mtx to long format for SummaryM.pool
# x) For SummaryM.pool how to define predict function?

# **********
# TO DO (ContinSummaryModel, binirize)
# * (***BUG***) Currently make.bins_mtx_1 fails on binary A with automatic bin detection (all values are placed in last (2nd) bin)
# * Need to handle gvars$misval values for self$outvar (contin) when transforming to bins (in fit, predict, predictAeqa)...
# * For x_cat the process is identical, except that normalize, define.intervals & make.ordinal is skipped
# * Need to test that make.bins_mtx_1 will do the right thing when x_cat is (0, ..., ncats) instead of (1, ..., ncats) (IT SHOULD WORK)

RegressionClass <- R6Class("RegressionClass",
	class = TRUE,
	portable = TRUE,
	public = list(
    useglm = TRUE,              # (Future vs.) TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
		family = NULL, 				      # (Future vs.) to run w/ other than "binomial" family	
		outvar.class = character(), # class of the outcome var: bin / cont / cat
		outvar = character(),       # regression outcome variable name
		predvars = character(),     # vector of regression covariate names (predictors)
		subset = NULL,              # subset expression (later evaluated to logical vector in the envir of the data)
		nbins = NULL,               # for cont. outvar, copied by ContinSummaryModel$new from O.datnetA$nbins.sVar(reg$outvar)
		bin_nms = NULL,
		pool.cont = FALSE,          # pool observations (long format) for binarized continous outcomes from different bins
		initialize = function(outvar.class = gvars$sVartypes$bin, outvar, predvars, subset, useglm = TRUE) {
			self$outvar.class <- outvar.class
			class(self) <- c(self$outvar.class, class(self))
			self$outvar <- outvar
			self$predvars <- predvars
			self$subset <- subset
			self$useglm <- useglm
		}
	),
	active = list(
		reg = function() {
			list(outvar.class = self$outvar.class,
				outvar = self$outvar,
				predvars = self$predvars,
				subset = self$subset)
		}
	)
)

NewSummaryModel <- function(reg, O.datnetA, ...) { UseMethod("NewSummaryModel") } # Generic S3 constructor for the summary model classes

NewSummaryModel.contin <- function(reg, O.datnetA, ...) { # Summary model constructor for continuous outcome sA[j]
	print("ContinSummaryModel constructor called...")
	ContinSummaryModel$new(reg = reg, O.datnetA = O.datnetA, ...)
	# Alternative name: ContOutModel
}

# NewSummaryModel.cat = function(reg, ...) { # Summary model constructor for categorical outcome sA[j]
#   print("CatSummaryModel constructor called...")
#   CatSummaryModel$new(...)
# }

# Summary model constructor for binary outcome sA[j]
NewSummaryModel.binary = function(reg, ...) {
	print("BinOutModel constructor called...")
	BinOutModel$new(glm = TRUE, reg = reg, ...) # fit a model with new object BinOutModel using glm.fit
	# BinOutModel$new(glm = FALSE, reg = reg, ...) # fit a model with new object BinOutModel using speedglm.wfit
	# Alternative name: BinOutModel
}


## ---------------------------------------------------------------------
# Class for defining, managing, fitting and returning the likelihood P(sA = sa | sW = sw) under g_star or g_0;
# Accepts (1) data (data.frame) for (sA,sW), (2) newdata (data.frame) for prediction, (3) obsdat.sA (matrix) for sa values;
# Defines and manages the factorization of the joint P(sA = sa | ... ) into reg models sA[j] ~ \bar{sA[j-1]} + sW;
# Figures out reg mdel factorization based on name ordering in (sA_nms, sW_nms);
# Evaluates subset_exprs in the envirs of data and newdata data.frames
# Calls BinOutModel$new, assumes each sA[j] is binary in reg (sA[j] ~ \bar{sA[j-1]} + sW);
## ---------------------------------------------------------------------

#' @title Class for defining, holding and fitting collections of summary measure models P(sA[j]|sW,\bar{sA}[j])
#' @docType class
#' @format An R6 class object.
#' @name SummariesModel
#' @details Following fields are created during initialization
#' \itemize{
#' \item{n_regs} ...
#' \item{nodes} ...
#' \item{sA_nms} ...
#' \item{sW_nms} ...
#' }
#' More details about the class...
##' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
##' @export
SummariesModel <- R6Class(classname = "SummariesModel",
	portable = TRUE,
	class = TRUE,
	public = list(
		n_regs = integer(),     # total no. of reg. models (logistic regressions)
		sA_nms = character(),   # sA names
		sW_nms = character(),   # sW names

		initialize = function(sA_class, sA_nms, sW_nms, subset, ...) {
			self$sA_nms <- sA_nms
			self$sW_nms <- sW_nms
			self$n_regs <- n_regs <- length(self$sA_nms)  # Number of sep. logistic regressions, determines the fitting algorithm this object will use (sep.logist by k, pooled, etc):
			print("#----------------------------------------------------------------------------------");
			print("New SummariesModel object:"); 
			print("No. of regressions: " %+% self$n_regs)
			print("sA_classes: "); str(sA_class)
			print("sA_nms: "); print(self$sA_nms)
			print("sW_nms: "); print(self$sW_nms)
			print("#----------------------------------------------------------------------------------");

			if (!missing(subset)) {
				if (length(subset) < n_regs) {
					subset <- rep_len(subset, n_regs)
				} else if (length(subset) > n_regs) {
					# ... TO FINISH ...
					# increase n_regs to all combinations of (n_regs x subset)
				}
			} else {
				subset <- rep_len(list(TRUE), n_regs)
			}

			for (k_i in (1:n_regs)) {  # factorization by nFriends (k) value
				sA_i_nm <- self$sA_nms[k_i] # A variable we are predicting
				covars_nms <- c(self$sA_nms[-c(k_i:n_regs)], self$sW_nms) # dependent covars

				# Changed reg to object of RegressionClass:
				reg <- RegressionClass$new(outvar.class = sA_class[[k_i]], outvar = sA_i_nm, predvars = covars_nms, subset = subset[[k_i]])
				PsAsW.model <- NewSummaryModel(reg = reg, ...) # Constructor for new summary model P(sA[j]|\bar{sA}[j-1], sW}) object
				private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
				names(private$PsAsW.models)[k_i] <- "P(sA|sW)."%+%k_i
			}
			invisible(self)
		},
		length = function(){ base::length(private$PsAsW.models) },
		getfitted.pbins = function() { private$fitted.pbins },  # get all summaries objs
		getcumprodAeqa = function() { private$cumprodAeqa },  # get a vector of cum prod of P(sA[j]=a[j]|sW)
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
		  print("length(private$PsAsW.models)"); print(length(private$PsAsW.models))
		  for (k_i in seq_along(private$PsAsW.models)) { # loop over all regressions in PsAsW.models
		    private$PsAsW.models[[k_i]]$predict(newdata = newdata)
		  }
		  invisible(self)
		},

		# WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
		# Use daughter objects (stored from prev call to fit()) to run predict on P(sA=obsdat.sA|sW)
		# Invisibly return cumm. prob P(sA=sa|sW=sw)
		predictAeqa = function(obs.DatNet.sWsA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
			assert_that(is.DatNet.sWsA(obs.DatNet.sWsA))
			n <- obs.DatNet.sWsA$nobs
			cumprodAeqa <- rep_len(1, n)
			for (k_i in seq_along(private$PsAsW.models)) { # loop over all regressions in PsAsW.models
				print("getting the likelihood for outvar: "); print(private$PsAsW.models[[k_i]]$reg$outvar)
				cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA)
			}
			private$cumprodAeqa <- cumprodAeqa
			invisible(cumprodAeqa)
		}
	),

	active = list(
		actplhold = function() {}
	),

	private = list(
		PsAsW.models = list(),
		fitted.pbins = list(),
		cumprodAeqa = NULL
	)
)

## ---------------------------------------------------------------------
# Inherits from SummariesModel, allowing sA[j] in sA to be a continuous/categorical summary measure
# Defines the fitting algorithm for sA[j] ~ \bar{sA[j-1]} + sW
# Reconstructs the likelihood P(sA[j]=sa[j] | sW)
# If needed discretizes the continuous sA[j] and correctly defines regressions BinsA[j][i] ~ \bar{BinsA[j][i-1]} + \bar{sA[j-1]} + sW
# Creates the appropriate dataset of discretized summary measures (BinsA[j][1],...,BinsA[j][M])?
# Called from SummariesModel for contin sA[j]. Gets passed new subset definitions, e.g., (!mis(BinsA))
## ---------------------------------------------------------------------
ContinSummaryModel <- R6Class(classname = "ContinSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    O.datnetA = NULL,
    reg = NULL,
    outvar = character(),     # the name of the continous outcome var (sA[j])

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
        O.datnetA$set.sVar.intrvls(reg$outvar, O.datnetA$detect.sVar.intrvls(reg$outvar))
        self$reg$nbins <- O.datnetA$nbins.sVar(reg$outvar)
      }

      self$reg$bin_nms <- O.datnetA$bin.nms.sVar(reg$outvar) # new.sA_nms <- O.datnetA$bin.nms.sVar(reg$outvar)
      new.sA_class <- as.list(rep_len(gvars$sVartypes$bin, self$reg$nbins))
      names(new.sA_class) <- self$reg$bin_nms
      bin.m.params <- list(sA_class = new.sA_class, sA_nms = self$reg$bin_nms, sW_nms = reg$predvars)

      # *) new.subsets:
      # #todo 27 (ContinSummaryModel) +0: Put subset eval in a separate function (with var as arg + additional args) +
      # move new.subsets def into another location (inside ContinSummaryModel$new()?)
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
      bin.m.params <- append(bin.m.params, list(subset = new.subsets))

      # Combine list of params binparams from O.datnetA with additinal params passed in ...:
      addl_params <- list(...)
      bin.m.params <- append(bin.m.params, addl_params)
      parnames <- names(bin.m.params)
      if (length(bin.m.params) != 0 && (is.null(parnames) || any(parnames==""))) {
        stop("need to specify the name of each attribute")
      }

      print("new in continuous summary...")
      print("contin sA name: "%+%self$outvar);
      print("contin sA reg$nbins"); print(self$reg$nbins)
      print("contin binned sA names: "); print(self$reg$bin_nms)
      print("Contin. new.subsets: "); str(new.subsets)

      do.call(super$initialize, bin.m.params)  # call the parent class contructor
      # super$initialize(sA_class = new.sA_class, sA_nms = self$reg$bin_nms, sW_nms = new.sW_nms, subset = new.subsets, ...)
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
      # IF bins were defined inside ContinSummaryModel$new(), they will be saved in O.datnetA and DatNet.sWsA will not see them
      # unless the intervals are copied from O.datnetA
      data$copy.cbin.intrvls()

      print("intervals after data$copy.cbin.intrvls() in ContinSummaryModel$fit: "); print(data$get.sVar.intrvls(self$outvar))
      print("data$dat.sVar for: " %+% self$outvar); print(head(data$dat.sVar))
      print("currently active binirized sVar: " %+% data$active.bin.sVar)

      # Evaluates subsets & saves correct matrix to DatNet.sWsA
      data$binirize.sVar(self$outvar) # Note, don't need to save mat_bin here, its already saved inside DatNet.sWsA
      # mat_bin <- data$binirize.sVar(self$outvar)
      print("binirized data for: " %+% self$outvar)
      print("currently active binirized sVar: " %+% data$active.bin.sVar)

      print(head(data$binirize.sVar(self$outvar)))

      super$fit(data) # call the parent class fit method
      message("fit for " %+% self$outvar %+% " var succeeded...")
      invisible(self)
    },

    predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
      if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
        return(invisible(self))
      }
      print("predict in continuous summary...")
      # NEW VERSION. mat_bin doesn't need to be saved (even though its invisibly returned). 
      # mat_bin is automatically saved in datnet.sW.sA - potentially dangerous!!!
      mat_bin <- newdata$binirize.sVar(self$outvar)
      # print("mat_bin"); print(head(mat_bin))
      super$predict(newdata)
      invisible(self)
    },

    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(obs.DatNet.sWsA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      print("predictAeqa in continuous summary...")

      assert_that(is.DatNet.sWsA(obs.DatNet.sWsA))
      assert_that(is.vector(obs.DatNet.sWsA$get.outvar(var = self$outvar)))
      obs.DatNet.sWsA$binirize.sVar(self$outvar) # mat_bin <- obs.DatNet.sWsA$binirize.sVar(self$outvar)
      cumprodAeqa <- super$predictAeqa(obs.DatNet.sWsA = obs.DatNet.sWsA)

      # cumprodAeqa <- super$predictAeqa(self$binirize(obsdat.sA)) # old way
      invisible(cumprodAeqa)
      # invisible(super$predictAeqa(self$binirize(obsdat.sA)))  # one line alternative to the above?
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)

## ---------------------------------------------------------------------
# Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] 
# add a covariate BinInd (1,..,M), convert to long format (turn BinsA[1], ...., BinsA[M] into rows) 
# then call $super$fit on transformed data?
## ---------------------------------------------------------------------
ContinSummaryModelPool <- R6Class(classname = "ContinSummaryModelPool",
  inherit = ContinSummaryModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
    initialize = function(...) {
    	super$initialize(...) # call the parent class contructor
    }
  ),
  active = list(
    actplhold = function() {}
  ),
  private = list(
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)

# THIS IS A BAD APPROACH (DEFINING SEPARATE CLASS FOR CAT). NEED TO HANDLE cat sA[j] INSIDE ContinSummaryModel
CatSummaryModel <- R6Class(classname = "CatSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE
)
