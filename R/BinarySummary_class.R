


logitlinkinv <- function (eta) .Call(stats:::"C_logit_linkinv", eta) # use this glm logitlink inverse function

logisfit = function(datsum_obj) { UseMethod("logisfit") } # Generic for fitting the logistic model

logisfit.glmS3 = function(datsum_obj) { # # S3 method for glm binomial family fit, takes DatBinarySummary data object 
  message("calling glm.fit...")
  Xmat = datsum_obj$getXmat 
  Y_vals = datsum_obj$getY
  if (nrow(Xmat) == 0L) { # Xmat has 0 rows: return NA's and avoid throwing exception
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

logisfit.speedglmS3 = function(datsum_obj) { # S3 method for speedglm binomial family fit, takes DatBinarySummary data object 
  message("calling speedglm.wfit...")
  Xmat = datsum_obj$getXmat
  Y_vals = datsum_obj$getY
  print("dims of glm data:")
  print(dim(Xmat)); print(length(Y_vals))
  if (nrow(Xmat) == 0L) { # Xmat has 0 rows: return NA's and avoid throwing exception
    m.fit <- list(coef = rep_len(NA_real_, ncol(Xmat)))
  } else {
    m.fit <- speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial())
  }
  fit <- list(coef = m.fit$coef, linkfun = logitlinkinv, fitfunname = "speedglm")
  class(fit) <- c(class(fit), c("speedglmS3"))
  return(fit)
}

## ---------------------------------------------------------------------
#' Abstract summary measure class for P(sA[j]|sW,\bar{sA}[j])
#'
#' @importFrom R6 R6Class
#' @export
Abstract_DatBinarySummary <- R6Class(classname = "Abstract_DatBinarySummary",
# Abstract_DatBin <- R6Class(classname = "Abstract_DatBin",  
  portable = TRUE,
  class = TRUE,
  public = list(
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    n = NA_integer_,        # number of rows in the input data
    subset_expr = NULL,     # PASS THE LOGICAL EXPRESSIONS TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)
    initialize = function(...) { stop("cannot create abstract data store directly")},
    show = function() { # printing regression
      "P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=",") %+% ")"
    }
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
#' Data storage class for one summary measure (sA[j]|sW,\bar{sA}[j]).
#' To be used for storing and passing the design mat + outcome(s) + whatever else is needed.
#' Could include a reference to the full data.frame (as a field)
#' Should be able to store/create data in various ways 
#' Should be able to convert wide to long
#' Should be able to subset
#' Will store methods that perform necessary data transformations
#' Consider making this just a closure (class = FALSE) for faster dispatch, in which case S3 dispatch on DatBinarySummary object will not be possible
#' Add method to save design matrix when its created? Second time getXmat(), return previously created dmat?
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.string is.flag
#' @export

# TO DO: rename class to DatBin
DatBinarySummary <- R6Class(classname = "DatBinarySummary",
# DatBin <- R6Class(classname = "DatBin",
  inherit = Abstract_DatBinarySummary,
  # inherit = Abstract_DatBin,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(reg, ...) {
      assert_that(is.string(reg$outvar))
      assert_that(is.character(reg$predvars))
      self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      self$subset_expr <- reg$subset
      if (is.null(self$subset_expr)) {self$subset_expr <- TRUE}
        # self$subset_idx <- rep_len(TRUE, nrow(data))
        # self$subset_idx <- as.expression(TRUE) # an alternative to reduce the no. of checks on self$subset_idx
        # self$subset_idx <- as.call(TRUE) # an alternative to reduce the no. of checks on self$subset_idx
      assert_that(is.logical(self$subset_expr) || is.call(self$subset_expr))
      # print("Declared new DatBinarySummary"); print(self$show()); print("with subset expr: "); print(self$subset_expr)
      # self$setdata(data = data, subset = subset, ...)
      invisible(self)
    },

    # TO DO: move to private method...
    # Use data$get.df.sW.sA(covars, subset) to get data.frame of all predictors, data$vec.var(var, subset) to get outcome var vector 
    # Change below to: setdata = function(dat.sWsA, ...) { # Sets X_mat, Yvals, evaluates subset and performs correct subseting of data
    setdata = function(data, ...) { # Sets X_mat, Yvals, evaluates subset and performs correct subseting of data
      assert_that(is.DatNet.sWsA(data))
      self$n <- data$nobs

      if (is.logical(self$subset_expr)) {self$subset_idx <- self$subset_expr}
      if (is.call(self$subset_expr)) {
        # THIS IS GOING TO BE DIFFICULT. 
        # evalsubst needs to be changed so that environments of both dat.sA and dat.sW are available
        self$subset_idx <- data$evalsubst(subsetexpr = self$subset_expr)
        # change to: self$subset_idx <- datsWsA$evalsubst(subsetexpr = self$subset_expr)
        # self$subset_idx <- evalsubst(data = data, subsetexpr = self$subset_expr)

        print("self$subset_expr: "); print(self$subset_expr)
        print("fit subset_idx: "%+%length(self$subset_idx)); 
        print("sum subset: "%+%sum(self$subset_idx)); print(head(self$subset_idx))

        assert_that(is.logical(self$subset_idx))
        assert_that((length(self$subset_idx) == self$n) || (length(self$subset_idx) == 1L))
        # e.g., subset <- TRUE means select all rows or subset <- "(nFriends==3)"
      }
      private$Y_vals <- data$get.outvar(self$subset_idx, self$outvar) # Will always be a vector

      if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
        private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
        colnames(private$X_mat) <- c("Intercept", self$predvars)
      } else {
        print("Dat Bin self$predvars"); print(self$predvars)
        private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.df.sW.sA(self$subset_idx, self$predvars, self$outvar)))
        # old version for X_mat when data was an actual data.frame of all cbind(sW,sA)
        # private$X_mat <- as.matrix(cbind(Intercept = 1, data[self$subset_idx, self$predvars, drop = FALSE]))
        # To find and replace misvals in X_mat: # private$X_mat[private$X_mat == gvars$misval] <- gvars$misXreplace
      }
      print("Y_vals: "); print(length(private$Y_vals)); print(head(private$Y_vals))
      print("design mat: "); print(dim(private$X_mat)); print(class(private$X_mat)); print(head(private$X_mat))
    },

    newdata = function(newdata, ...) {
      assert_that(is.DatNet.sWsA(newdata)) # old: assert_that(is.data.frame(newdata)||is.matrix(newdata))
      self$setdata(data = newdata, ...)
      invisible(self)
    },

    # Generic prediction fun for logistic regression coefs, predicts P(A=1|newXmat)
    # No need for S3 for now, until need different pred. funs for different classes
    # Does not handle cases with deterministic Anodes in the original data..
    logispredict = function(m.fit) {
      assert_that(!is.null(private$X_mat)); assert_that(!is.null(self$subset_idx))
      pAout <- rep_len(gvars$misval, self$n) # Set to default missing value for A[i] degenerate/degerministic/misval: # Alternative, set to default replacement val: pAout <- rep_len(gvars$misXreplace, newdatsum_obj$n)      
      if (sum(self$subset_idx > 0)) {
        eta <- private$X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
        pAout[self$subset_idx] <- m.fit$linkfun(eta)
      }
      return(pAout)
    }
    ),

  active = list( # 2 types of active bindings (w and wout args)
   emptydata = function() {
      self$n <- NA_integer_
      self$subset_idx <- NULL # wipe out subset as well or leave it alone?
      private$Y_vals <- NULL
      private$X_mat <- NULL
      
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
#' Class for modeling one summary measure, P(sA[j]|sW,\bar{sA}[j])
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
#' @export

BinarySummaryModel <- R6Class(classname = "BinarySummaryModel",
# BinOutModel  <- R6Class(classname = "BinOutModel",
  inherit = DatBinarySummary,
  # inherit = DatBin,
  portable = TRUE,
  class = TRUE,
  public = list(
    glmfitclass = "glmS3", # default glm fit class
    is.fitted = FALSE,
    datbin = NULL, # object of class DatBinarySummary that is used in fitting / prediction, never saved (need to be initialized with $new())

    initialize = function(glm = TRUE, reg, ...) {
      assert_that(is.flag(glm)) # THIS MIGHT BE MOVED TO self$fit() FUNCTION
      if (!glm) self$glmfitclass <- "speedglmS3"
      self$datbin <- DatBinarySummary$new(reg = reg, ...) # postponed adding data in DatBinarySummary until self$fit() is called
      class(self$datbin) <- c(class(self$datbin), self$glmfitclass)
      # can also use: self$datbin <- DatBinarySummary$new(glm = self$glmfitclass, ...)
      # or: self$datbin <- DatBinarySummary$new(self, ...) (passing self might get confusing)
      print("New BinarySummaryModel:"); print(self$show())
      invisible(self)
    },

    show = function() {self$datbin$show()},
    getfit = function() {private$m.fit},
    getprobA1 = function() {private$probA1},
    getprobAeqa = function() {private$probAeqa},

    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      # replace with: if (!self$overwrite) assert_that(!self$is.fitted)
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      self$datbin$newdata(newdata = data, ...) # populate datbin with X_mat & Y_vals
      # ... additional checks & assertions... # ... add try() to below: # ... add checks for testing a successful fit

      private$m.fit <- logisfit(datsum_obj = self$datbin) # private$m.fit <- data_obj$logisfit or private$m.fit <- data_obj$logisfit() # alternative 2 is to apply data_obj method / method that fits the model
      # print("private$m.fit"); print(private$m.fit)
      self$is.fitted <- TRUE

      private$probA1 <- self$datbin$logispredict(m.fit = private$m.fit)
      # private$probA1 <- logispredict(newdatsum_obj = self$datbin, m.fit = private$m.fit)
      # print("fit private$m.fit"); print(private$m.fit)
      # print("fit private$probA1"); print(private$probA1)

      self$datbin$emptydata  # Xmat, Yvals & subset in datbin no longer needed
      invisible(self)
    },

    predict = function(newdata, ...) { # P(A^s[i]=1|W^s=w^s): uses private$m.fit to generate predictions for newdata
      assert_that(self$is.fitted) # vs. stopifnot(self$is.fitted)
      if (missing(newdata)) { # ... Do nothing, predictions for fitted data are already saved
        return(invisible(self))
      }
      self$datbin$newdata(newdata = newdata, ...) # re-populate datbin with new X_mat & Y_vals
      # print("predict newdata"); print(head(newdata))
      private$probA1 <- self$datbin$logispredict(m.fit = private$m.fit) # overwrite probA1 with new predictions:
      # private$probA1 <- logispredict(newdatsum_obj = self$datbin, m.fit = private$m.fit) # overwrite probA1 with new predictions:
      # print("predict private$probA1: "); print(private$probA1)
      invisible(self) # returning self allows chaining object manipulations
    },
    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(obsdat.sA) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a's)
      assert_that(self$is.fitted)
      assert_that(is.integer(obsdat.sA)) # check B: obsdat.sA is a row of integers
      assert_that(is.logical(self$getsubset))

      # print("private$probA1"); print(private$probA1)
      # print("sum(self$getsubset)"); print(sum(self$getsubset)); print(self$getsubset)
      # print("private$probA1[self$getsubset]"); print(private$probA1[self$getsubset])

      assert_that(!any(is.na(private$probA1[self$getsubset]))) # check that predictions P(A=1|dmat) exist for all obs.

      private$probAeqa <- rep_len(1L, length(obsdat.sA)) # for missing, the likelihood is always set to P(A = a) = 1.
      # TO FIND AND REPLACE FOR MISSING VALS IN obsdat.sA: # obsdat.sA[obsdat.sA==gvars$misval] <- gvars$misXreplace

      probA1 <- private$probA1[self$getsubset];
      indA <- obsdat.sA[self$getsubset]
      private$probAeqa[self$getsubset] <- probA1^(indA) * (1L - probA1)^(1L - indA)
      # print("private$probA1: "); print(private$probA1)
      # print("self$getsubset: "); print(self$getsubset)
      # print("obsdat.sA: "); print(obsdat.sA)
      # print("private$probAeqa: "); print(private$probAeqa)
      self$datbin$emptydata  # wipe out prediction data after getting the likelihood
      invisible(private$probAeqa)
    }
  ),

  active = list(
    getsubset = function() {self$datbin$subset_idx}
    # x2 = function(value) {
    #   if (missing(value)) return(self$x * 2)
    #   else self$x <- value/2
    # },
    # rand = function() {rnorm(1)}
  ),

  private = list(
    m.fit = list(), # the model fit (coefficients)
    probA1 = NULL,    # Predicted probA^s=1 conditional on X_mat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on X_mat
  )
)
