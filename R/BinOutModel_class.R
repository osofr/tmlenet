#-----------------------------------------------------------------------------
# TO DO: 
#-----------------------------------------------------------------------------
# - #todo 35 (NewSummaryModel.binary) +0: Change constructor arg glm to be always inside reg object
# - (Low priority) Consider merging these two classes into one (BinDat & BinOutModel)

logitlinkinv <- function (eta) .Call(stats:::"C_logit_linkinv", eta) # glm logitlink inverse fun
logisfit <- function(datsum_obj) UseMethod("logisfit") # Generic for fitting the logistic model
# S3 method for glm binomial family fit, takes BinDat data object:
logisfit.glmS3 <- function(datsum_obj) {
	message("calling glm.fit...")
	Xmat <- datsum_obj$getXmat
	Y_vals <- datsum_obj$getY
    # Xmat has 0 rows: return NA's and avoid throwing exception:
	if (nrow(Xmat) == 0L) {
	  m.fit <- list(coef = rep_len(NA_real_, ncol(Xmat)))
	} else {
	  ctrl <- glm.control(trace=FALSE, maxit=1000)
	  # SuppressGivenWarnings({
	            m.fit <- glm.fit(x = Xmat, y = Y_vals, family = binomial(), control = ctrl)
	            # }, GetWarningsToSuppress())
	}
	fit <- list(coef = m.fit$coef, linkfun = logitlinkinv, fitfunname = "glm")
	class(fit) <- c(class(fit), c("glmS3"))
	return(fit)
}

logisfit.speedglmS3 <- function(datsum_obj) { # S3 method for speedglm binomial family fit, takes BinDat data object
	# message("calling speedglm.wfit...")
	Xmat <- datsum_obj$getXmat
	Y_vals <- datsum_obj$getY

	# print("dims of glm data:")
	# print(dim(Xmat));
	# print("head(Xmat)");; print(head(Xmat));
	# print(length(Y_vals));
	# print(head(Y_vals))

	if (nrow(Xmat) == 0L) { # Xmat has 0 rows: return NA's and avoid throwing exception
	  m.fit <- list(coef = rep_len(NA_real_, ncol(Xmat)))
	} else {
	  m.fit <- speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial())
	  # m.fit <- speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial(), sparse = TRUE)
	}
	fit <- list(coef = m.fit$coef, linkfun = logitlinkinv, fitfunname = "speedglm")
	# print("fit"); print(fit)

	class(fit) <- c(class(fit), c("speedglmS3"))
	return(fit)
}

# S3 methods for getting coefs from fitted BinOutModel class object
coef.BinOutModel <- function(binoutmodel) {
	assert_that(binoutmodel$is.fitted)
	fit <- binoutmodel$getfit
	fit$coef
}
summary.BinOutModel <- function(binoutmodel) {
	assert_that(binoutmodel$is.fitted)
	fit <- binoutmodel$getfit
	append(list(reg = binoutmodel$show()), fit)
}

#-----------------------------------------------------------------------------
# TO DO: Incorporate this inside BinDat / BinOutModel classes, to be called when reg object field !is.null(form)
# Formula based glm model fit 
#-----------------------------------------------------------------------------
f_est <- function(d, form, family) {
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
    SuppressGivenWarnings({
              m <- glm(as.formula(form),
                  data = d,
                  family = family,
                  control = ctrl)
              },
              GetWarningsToSuppress())
    return(m)
}

## ---------------------------------------------------------------------
#' (NOT USED) Abstract summary measure class for P(sA[j]|sW,sA[j]) 
#'
#' @export
Abstract_BinDat <- R6Class(classname = "Abstract_BinDat",
	portable = TRUE,
	class = TRUE,
	public = list(
	  initialize = function(...) { stop("cannot create abstract data store directly")}
	  # store = function(..., key = digest(list(...)), overwrite = FALSE,
	  #   envir = parent.frame()) {
	  #   stop("store method not implemented")
	  # },
	  # try_load = function(key, envir = parent.frame()) {
	  #   ads_try_load(self, key, envir)
	  # },
	)
)
# ## ---------------------------------------------------------------------
# USAGE:
# #' @keywords internal
# ads_try_load <- function(self, key, envir) {
#   try(self$load(key, envir = envir), silent = TRUE)
# }

## ---------------------------------------------------------------------
#' Data storage class for one summary measure (sA[j] , sW,sA[j]).
#' To be used for storing and passing the design mat + outcome(s) + whatever else is needed.
#' Could include a reference to the full data.frame (as a field)
#' Should be able to store/create data in various ways 
#' Should be able to convert wide to long
#' Should be able to subset
#' Will store methods that perform necessary data transformations
#' Consider making this just a closure (class = FALSE) for faster dispatch, in which case S3 dispatch on BinDat object will not be possible
#' Add method to save design matrix when its created? Second time getXmat(), return previously created dmat?
#'
#' @importFrom assertthat assert_that is.count is.string is.flag
#' @export

BinDat <- R6Class(classname = "BinDat",
	# inherit = Abstract_BinDat,
	portable = TRUE,
	class = TRUE,
	public = list(
		reg = NULL,
		outvar = character(),   # outcome name(s)
		predvars = character(), # names of predictor vars
		n = NA_integer_,        # number of rows in the input data
		subset_expr = NULL,     # PASS THE LOGICAL EXPRESSIONS TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
		subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)
		initialize = function(reg, ...) {
			assert_that(is.string(reg$outvar))
			assert_that(is.character(reg$predvars))
			self$reg <- reg
			self$outvar <- reg$outvar
			self$predvars <- reg$predvars
			self$subset_expr <- reg$subset
			if (is.null(self$subset_expr)) {self$subset_expr <- TRUE}
			assert_that(is.logical(self$subset_expr) || is.call(self$subset_expr))
			# print("Declared new BinDat"); print(self$show()); print("with subset expr: "); print(self$subset_expr)
			invisible(self)
		},
        # printing regression:
		show = function() {
			"P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=",") %+% ")"
		},

		newdata = function(newdata, getoutcome = FALSE, ...) {
			assert_that(is.DatNet.sWsA(newdata)) # old: assert_that(is.data.frame(newdata)||is.matrix(newdata))
			self$setdata(data = newdata, getoutcome = getoutcome, ...)
			invisible(self)
		},

		# TO DO: move to private method...
        # Sets X_mat, Yvals, evaluates subset and performs correct subseting of data
		setdata = function(data, getoutcome, ...) {
			assert_that(is.DatNet.sWsA(data))
			self$n <- data$nobs
			if (is.logical(self$subset_expr)) {self$subset_idx <- self$subset_expr}
			if (is.call(self$subset_expr)) {
				self$subset_idx <- data$evalsubst(subsetexpr = self$subset_expr)
				assert_that(is.logical(self$subset_idx))
				assert_that((length(self$subset_idx) == self$n) || (length(self$subset_idx) == 1L))
				# e.g., subset <- TRUE means select all rows or subset <- "(nFriends==3)"
			}

			if (getoutcome) private$Y_vals <- data$get.outvar(self$subset_idx, self$outvar) # Always a vector

			if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
				private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
				colnames(private$X_mat) <- c("Intercept", self$predvars)
			} else {
				private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.dat.sWsA(self$subset_idx, self$predvars)))
				# To find and replace misvals in X_mat:
				if (self$reg$ReplMisVal0) {
					private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace
				}
			}
		},

		# Generic prediction fun for logistic regression coefs, predicts P(A=1|newXmat)
		# No need for S3 for now, until need different pred. funs for different classes
		# Does not handle cases with deterministic Anodes in the original data..
		logispredict = function(m.fit) {
			assert_that(!is.null(private$X_mat)); assert_that(!is.null(self$subset_idx))
            # Set to default missing value for A[i] degenerate/degerministic/misval:
            # Alternative, set to default replacement val: pAout <- rep_len(gvars$misXreplace, newdatsum_obj$n)

			pAout <- rep_len(gvars$misval, self$n)
			# 07/14/15: Need probA1 even for degenerate bins to be able to normalize by bin-width correctly
			# pAout <- rep_len(1L, self$n)

			if (sum(self$subset_idx > 0)) {
				eta <- private$X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
				pAout[self$subset_idx] <- m.fit$linkfun(eta)
			}
			return(pAout)
		}
	),

	active = list( # 2 types of active bindings (w and wout args)
		emptydata = function() {
			private$X_mat <- NULL
			# self$n <- NA_integer_
			# private$Y_vals <- NULL # Can't wipe out the subset and make predictions for missing(newdata)
			# self$subset_idx <- NULL # Can't wipe out the subset and make predictions for missing(newdata)
		},
		getXmat = function() {private$X_mat},
		getY = function() {private$Y_vals}
		),

	private = list(
		X_mat = NULL,
		Y_vals = NULL
	 )
	)

## ---------------------------------------------------------------------
#' Class for fitting a logistic model with binary outcome, P(sA[j]|sW,sA[j])
#'
#' @importFrom assertthat assert_that is.flag
#' @export

BinOutModel  <- R6Class(classname = "BinOutModel",
	portable = TRUE,
	class = TRUE,
	public = list(
		reg = NULL,
		glmfitclass = "glmS3", # default glm fit class
		is.fitted = FALSE,
		bindat = NULL, # object of class BinDat that is used in fitting / prediction, never saved (need to be initialized with $new())

		# TO DO: Change to glm arg being inside reg (reg$glm):
		initialize = function(glm = TRUE, reg, ...) {
			assert_that(is.flag(glm)) # THIS MIGHT BE MOVED TO self$fit() FUNCTION
			if (!glm) self$glmfitclass <- "speedglmS3"
			self$reg <- reg
			self$bindat <- BinDat$new(reg = reg, ...) # postponed adding data in BinDat until self$fit() is called
			class(self$bindat) <- c(class(self$bindat), self$glmfitclass)
			# can also use: self$bindat <- BinDat$new(glm = self$glmfitclass, ...)
			# or: self$bindat <- BinDat$new(self, ...) (passing self might get confusing)
			print("Init BinOutModel:"); print(self$show())
			invisible(self)
		},

		fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
			if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked
			self$bindat$newdata(newdata = data, getoutcome = TRUE, ...) # populate bindat with X_mat & Y_vals
			# ... additional checks & assertions... # ... add try() to below: # ... add checks for testing a successful fit
			private$m.fit <- logisfit(datsum_obj = self$bindat) # private$m.fit <- data_obj$logisfit or private$m.fit <- data_obj$logisfit() # alternative 2 is to apply data_obj method / method that fits the model
			self$is.fitted <- TRUE
			private$probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
			self$bindat$emptydata  # Xmat in bindat is no longer needed, only subset, outvar & probA1 (probOut1)
			invisible(self)
		},

		# #todo 13 (BinOutModel, predict, predictAeqa) +0: Need to be linked together, since can create a discrepancy for missing(newdata) but !missing(obs.DatNet.sWsA)
		predict = function(newdata, ...) { # P(A^s[i]=1|W^s=w^s): uses private$m.fit to generate predictions for newdata
			assert_that(self$is.fitted) # vs. stopifnot(self$is.fitted)
            # ... Do nothing, predictions for fitted data are already saved
			if (missing(newdata)) {
				# stop("must provide newdata for BinOutModel$predict()")
				return(invisible(self))
			}
			self$bindat$newdata(newdata = newdata, getoutcome = FALSE, ...) # re-populate bindat with new X_mat & Y_vals
			private$probA1 <- self$bindat$logispredict(m.fit = private$m.fit) # overwrite probA1 with new predictions:
			invisible(self) # returning self allows chaining object manipulations
		},

		# WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
		# Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
		# Invisibly return cumm. prob P(sA=sa|sW=sw)

		# #todo 15 (BinOutModel, predictAeqa) +0: TO DO: MOVE PART OF THE CODE TO self
		predictAeqa = function(obs.DatNet.sWsA) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a's)
			assert_that(self$is.fitted)
			assert_that(is.logical(self$getsubset))

			# obtain predictions (likelihood) for response on fitted data:
			if (missing(obs.DatNet.sWsA)) {
				n <- length(self$getsubset)
				indA <- self$getoutvarval[self$getsubset] # Always a vector of 0/1
			} else {
				n <- obs.DatNet.sWsA$nobs
				indA <- obs.DatNet.sWsA$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
			}

			private$probAeqa <- rep_len(1L, n) # for missing, the likelihood is always set to P(A = a) = 1.
			assert_that(!any(is.na(private$probA1[self$getsubset]))) # check that predictions P(A=1|dmat) exist for all obs.
			probA1 <- private$probA1[self$getsubset];
			assert_that(is.integerish(indA)) # check B: obsdat.sA is a row of integers
			private$probAeqa[self$getsubset] <- probA1^(indA) * (1L - probA1)^(1L - indA)
			# private$probAeqa[self$getsubset] <- probA1^(indA)
			self$bindat$emptydata  # wipe out prediction data after getting the likelihood
			invisible(private$probAeqa)

		},
		
		show = function() {self$bindat$show()},
		setprobA1 = function(probA1) { private$probA1 <- probA1 }
	),
	active = list(
		getfit = function() { private$m.fit },
		getprobA1 = function() { private$probA1 },
		# getprobAeqa = function() { private$probAeqa },
		getsubset = function() { self$bindat$subset_idx },
		getoutvarval = function() { self$bindat$getY },
		getoutvarnm = function() { self$bindat$outvar }
	),
	private = list(
		m.fit = list(), # the model fit (coefficients)
		probA1 = NULL,    # Predicted probA^s=1 conditional on X_mat
		probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on X_mat
	)
)