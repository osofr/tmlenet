#-----------------------------------------------------------------------------
# TO DO: 
# - #todo 81 (DatBin, BitOutModel) +0: Create a new class DatBinLong that inherits from DatBin
 # with methods for setdata, predict overridden with pooled fit based versions
# - (Low priority) Consider merging these two classes into one (BinDat & BinOutModel)
#-----------------------------------------------------------------------------

logitlinkinv <- function (eta) .Call(stats:::"C_logit_linkinv", eta) # glm logitlink inverse fun
logisfit <- function(datsum_obj) UseMethod("logisfit") # Generic for fitting the logistic model

# S3 method for glm binomial family fit, takes BinDat data object:
logisfit.glmS3 <- function(datsum_obj) {
  if (gvars$verbose) message("calling glm.fit...")
  Xmat <- datsum_obj$getXmat
  Y_vals <- datsum_obj$getY
    # Xmat has 0 rows: return NA's and avoid throwing exception:
  if (nrow(Xmat) == 0L) {
    m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    ctrl <- glm.control(trace=FALSE, maxit=1000)
    # SuppressGivenWarnings({
              m.fit <- stats::glm.fit(x = Xmat, y = Y_vals, family = binomial(), control = ctrl)
              # }, GetWarningsToSuppress())
  }
  fit <- list(coef = m.fit$coef, linkfun = "logitlinkinv", fitfunname = "glm")

  if (gvars$verbose) print(fit$coef)

  class(fit) <- c(class(fit), c("glmS3"))
  return(fit)
}

# S3 method for speedglm binomial family fit, takes BinDat data object:
logisfit.speedglmS3 <- function(datsum_obj) {
  if (gvars$verbose) message("calling speedglm.wfit...")
  Xmat <- datsum_obj$getXmat
  Y_vals <- datsum_obj$getY

  if (nrow(Xmat) == 0L) { # Xmat has 0 rows: return NA's and avoid throwing exception
    m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    m.fit <- try(speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial()), silent=TRUE)
    # m.fit <- speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial(), sparse = TRUE)
    if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
      message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", m.fit)
      return(logisfit.glmS3(datsum_obj))
    }
  }
  fit <- list(coef = m.fit$coef, linkfun = "logitlinkinv", fitfunname = "speedglm")
  if (gvars$verbose) print(fit$coef)
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


# Convert existing Bin matrix (Bin indicators) for continuous self$outvar into long format data.table with 3 columns:
# ID - row number; sVar_allB.j - bin indicators collapsed into one col; bin_ID - bin number identify for prev. columns
# automatically removed all missing (degenerate) bin indicators
binirized.to.DTlong = function(BinsDat_wide, binID_seq, ID, 
                                bin_names, pooled_bin_name, name.sVar) {
  # Convert Bin matrix into a data.table (without data.frame as intermediate), with new row ID column:
  DT_BinsDat_wide <- as.data.table(BinsDat_wide)[, c("ID") := ID, with = FALSE]
  setcolorder(DT_BinsDat_wide, c("ID", names(DT_BinsDat_wide)[-ncol(DT_BinsDat_wide)]))
  # melt into long format:
  sVar_melt_DT <- data.table:::melt.data.table(DT_BinsDat_wide,
                                                id.vars = "ID",
                                                measure.vars = bin_names,
                                                value.name = pooled_bin_name,
                                                variable.name = name.sVar,
                                                variable.factor = FALSE,
                                                na.rm = FALSE)

  print("sVar_melt_DT init: "); print(sVar_melt_DT)
  nbin_rep <- rep(binID_seq, each = nrow(BinsDat_wide))
  print("pooled_bin_name: " %+% pooled_bin_name);
  # 1) Add bin_ID; 2) remove a column with Bin names; 3) remove all rows with NA value for outcome (degenerate bins)
  sVar_melt_DT <- sVar_melt_DT[, c("bin_ID") := list(nbin_rep)][, name.sVar := NULL, with = FALSE][!is.na(get(pooled_bin_name))]
  data.table::setkeyv(sVar_melt_DT, c("ID", "bin_ID"))  # sort by ID, bin_ID to prepare for merge with predictors (sW)
  return(sVar_melt_DT)
}
# Prepare predictors (sW/X_mat) as data.table, adding row IDs for a join
# Join with sVar_melt_DT that is already in long format
# Need to check that created IDs match exactly for both datasets
join.Xmat = function(X_mat, sVar_melt_DT, ID) {
  nIDs <- length(unique(sVar_melt_DT[["ID"]]))
  assert_that(nIDs == nrow(X_mat))
  X_mat_DT <- data.table::as.data.table(X_mat)[, c("ID") := ID, with = FALSE]
  data.table::setkeyv(X_mat_DT, c("ID")) # sort by ID
  print("X_mat_DT: "); print(X_mat_DT)
  sVar_melt_DT <- sVar_melt_DT[X_mat_DT] # Merge long format (self$pooled_bin_name, binIDs) with predictors (sW)
  print("sVar_melt_DT[1:10,]"); print(sVar_melt_DT[1:10,])
  return(sVar_melt_DT)
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
#' Data storage / manipulation class for modeling summary measures P(sA[j] | sW,bar{sA}[j-1]).
#' Converts data wide to long for pooling across binary indicators (when fitting one model for several indicators)
#' Queries DatNet.sWsA to get appropriate data columns &  row subsets
#'
#' @importFrom assertthat assert_that is.count is.string is.flag
#' @export
BinDat <- R6Class(classname = "BinDat",
  # inherit = Abstract_BinDat,
  cloneable = FALSE,
  portable = TRUE,
  class = TRUE,
  public = list(
    # reg = NULL,
    bin_names = NULL,
    ID = NULL,
    pooled_bin_name = NULL,
    # binID_seq = NULL,
    nbins = integer(),
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    pool_cont = logical(),
    outvars_to_pool = character(),
    ReplMisVal0 = logical(),
    n = NA_integer_,        # number of rows in the input data
    subset_expr = NULL,     # PASS THE LOGICAL EXPRESSIONS TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)

    initialize = function(reg, ...) {
      assert_that(is.string(reg$outvar))
      assert_that(is.character(reg$predvars))
      # self$reg <- reg
      self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      self$subset_expr <- reg$subset
      self$pool_cont <- reg$pool_cont
      self$outvars_to_pool <- reg$outvars_to_pool
      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- self$reg$nbins
      if (is.null(self$subset_expr)) {self$subset_expr <- TRUE}
      assert_that(is.logical(self$subset_expr) || is.call(self$subset_expr) || is.character(self$subset_expr))
      # print("Declared new BinDat"); print(self$show()); print("with subset expr: "); print(self$subset_expr)
      invisible(self)
    },

    # printing regression:
    show = function() {
      "P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=",") %+% ")"
    },

    newdata = function(newdata, ...) {
      assert_that(is.DatNet.sWsA(newdata)) # old: assert_that(is.data.frame(newdata)||is.matrix(newdata))
      # CALL self$setdata.long() when: 1) self$pool_cont is TRUE & 2) more than one outvars_to_pool
      if (self$pool_cont && length(self$outvars_to_pool)>1) {
        self$setdata.long(data = newdata, ...)
      } else {
        self$setdata(data = newdata, ...)
      }
      invisible(self)
    },

    define.subset_idx = function(data) {
      if (is.logical(self$subset_expr)) {
        subset_idx <- self$subset_expr
      } else if (is.call(self$subset_expr)) {
        subset_idx <- data$evalsubst(subsetexpr = self$subset_expr)
      } else if (is.character(self$subset_expr)) {
        subset_idx <- data$evalsubst(subsetvars = self$subset_expr)
      }
      assert_that(is.logical(subset_idx))
      if ((length(subset_idx) < self$n) && (length(subset_idx) > 1L)) {
        message("subset_idx has smaller length than self$n; repeating subset_idx p times, for p: " %+% data$p)
        subset_idx <- rep.int(subset_idx, data$p)
        if (length(subset_idx) != self$n) stop("BinDat$define.subset_idx: self$n is not equal to nobs*p!")
      }
      assert_that((length(subset_idx) == self$n) || (length(subset_idx) == 1L))
      return(subset_idx)
    },

    # TO DO: move to private method...
    # Sets X_mat, Yvals, evaluates subset and performs correct subseting of data
    # everything is performed using data$ methods (data is of class DatNet.sWsA)
    setdata = function(data, ...) {
      assert_that(is.DatNet.sWsA(data))
      self$n <- data$nobs
      self$subset_idx <- self$define.subset_idx(data)

      private$Y_vals <- data$get.outvar(self$subset_idx, self$outvar) # Always a vector
      # if (getoutcome) private$Y_vals <- data$get.outvar(self$subset_idx, self$outvar) # Always a vector
      if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
        private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
        colnames(private$X_mat) <- c("Intercept", self$predvars)
      } else {
        # *** THIS IS THE ONLY LOCATION IN THE PACKAGE WHERE CALL TO DatNet.sWsA$get.dat.sWsA() IS MADE ***
        private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.dat.sWsA(self$subset_idx, self$predvars)))        
        # To find and replace misvals in X_mat:
        if (self$ReplMisVal0) private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace
      }      
      # print("self$predvars"); print(self$predvars)
      # print("head(private$X_mat)"); print(head(private$X_mat))
      invisible(self)
    },

    setdata.long = function(data, ...) {
      assert_that(is.DatNet.sWsA(data))
      self$n <- data$nobs
      self$subset_idx <- self$define.subset_idx(data)
      if (!(data$active.bin.sVar %in% self$outvar)) { stop("currently binirized sVar does not match self$outvar argument") }

      # Setting up object fields related to pooling of continuous sA:
      self$pooled_bin_name <- data$pooled.bin.nm.sVar(self$outvar)
      self$bin_names <- self$outvars_to_pool
      binID_seq <- 1L:self$nbins
      BinsDat_wide <- data$get.dat.sWsA(self$subset_idx, self$outvars_to_pool)
      self$ID <- as.integer(1:nrow(BinsDat_wide))

      print("self$pooled_bin_name: " %+% self$pooled_bin_name)
      print("self$bin_names: "); print(self$bin_names)
      print("self$nbins: " %+% self$nbins);
      print("binID_seq: "); print(binID_seq)

      # To grab bin Ind mat directly (prob a bit faster): BinsDat_wide <- data$dat.bin.sVar[self$subset_idx, ]
      print("self$subset_idx: "); print(head(self$subset_idx))
      print("class(BinsDat_wide): " %+% class(BinsDat_wide));
      print("head(BinsDat_wide) in BinDat: "); print(head(BinsDat_wide))

      BinsDat_long <- binirized.to.DTlong(BinsDat_wide = BinsDat_wide, binID_seq = binID_seq, ID = self$ID, 
                                          bin_names = self$bin_names, pooled_bin_name = self$pooled_bin_name, 
                                          name.sVar = self$outvar)

      print("BinsDat_long: "); print(BinsDat_long)
      print("class(BinsDat_long): "); print(class(BinsDat_long))
      
      sVar_melt_DT <- join.Xmat(X_mat = data$get.dat.sWsA(self$subset_idx, self$predvars), 
                                sVar_melt_DT = BinsDat_long, ID = self$ID)

      # prepare design matrix for modeling w/ glm.fit or speedglm.wfit:
      X_mat <- sVar_melt_DT[,c("bin_ID", self$predvars), with=FALSE][, c("Intercept") := 1] # select bin_ID + predictors, add intercept column
      setcolorder(X_mat, c("Intercept", "bin_ID", self$predvars)) # re-order columns by reference (no copy)
      self$ID <- sVar_melt_DT[["ID"]]
      private$X_mat <- as.matrix(X_mat)
      print("private$X_mat[1:10,]"); print(private$X_mat[1:10,])
      private$Y_vals <- sVar_melt_DT[, self$pooled_bin_name, with = FALSE][[1]] # outcome vector:
      print("head(private$Y_vals)"); print(head(private$Y_vals, 100))

      # **************************************
      # TO FINISH...
      # **************************************
      # if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
      #   private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
      #   colnames(private$X_mat) <- c("Intercept", self$predvars)
      # } else {
      #   # *** THIS IS THE ONLY LOCATION IN THE PACKAGE WHERE CALL TO DatNet.sWsA$get.dat.sWsA() IS MADE ***
      #   private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.dat.sWsA(self$subset_idx, self$predvars)))
        # To find and replace misvals in X_mat:
        if (self$ReplMisVal0) private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace
      # }
    },

    logispredict.long = function(m.fit) {
      assert_that(!is.null(private$X_mat)); assert_that(!is.null(self$subset_idx))
      assert_that(nrow(private$X_mat)==length(private$Y_vals))
      pAout <- rep.int(gvars$misval, self$n)
      if (sum(self$subset_idx > 0)) {
        # -----------------------------------------------------------------
        # OBTAINING PREDICTIONS FOR LONG FORMAT P(Ind_j = 1 | Bin_j, W) BASED ON EXISTING POOLED FIT:
        # -----------------------------------------------------------------
        eta <- private$X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
        probA1 <- match.fun(FUN = m.fit$linkfun)(eta)
        # -----------------------------------------------------------------
        # GETTING ID-BASED PREDICTIONS (n) as cumprod of P(Ind_j = 1 | Bin_j, W) for j = 1, ..., K
        # -----------------------------------------------------------------
        ProbAeqa_long <- as.vector(probA1^(private$Y_vals) * (1L - probA1)^(1L - private$Y_vals))
        res_DT <- data.table(ID = self$ID, ProbAeqa_long = ProbAeqa_long)
        print("res_DT: "); print(res_DT)
        res_DT <- res_DT[, list(cumprob = cumprod(ProbAeqa_long)), by = ID]
        # print("res_DT w/ cumprob: "); print(res_DT)
        data.table::setkeyv(res_DT, c("ID")) # sort by ID
        res_DT_short <- res_DT[unique(res_DT[, key(res_DT), with = FALSE]), mult = 'last']
        print("res_DT_short: "); print(res_DT_short)
        ProbAeqa <- res_DT_short[["cumprob"]]

        print("length(ProbAeqa): " %+% length(ProbAeqa))
        print("head(ProbAeqa, 50)"); print(head(ProbAeqa, 50))
        pAout[self$subset_idx] <- ProbAeqa
      }
      return(pAout)
    },

    # Generic prediction fun for logistic regression coefs, predicts P(A = 1 | newXmat)
    # No need for S3 for now, until need different pred. funs for different classes
    # Does not handle cases with deterministic Anodes in the original data..
    logispredict = function(m.fit) {
      assert_that(!is.null(private$X_mat)); assert_that(!is.null(self$subset_idx))
      # Set to default missing value for A[i] degenerate/degerministic/misval:
      # Alternative, set to default replacement val: pAout <- rep.int(gvars$misXreplace, newdatsum_obj$n)
      pAout <- rep.int(gvars$misval, self$n)
      if (sum(self$subset_idx > 0)) {
        eta <- private$X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
        pAout[self$subset_idx] <- match.fun(FUN = m.fit$linkfun)(eta)
        # pAout[self$subset_idx] <- m.fit$linkfun(eta)
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
    emptyY = function() { private$Y_vals <- NULL},
    emptySubset_idx = function() { self$subset_idx <- NULL },
    emptyN = function() { self$n <- NA_integer_ },
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
  cloneable = FALSE,
  portable = TRUE,
  class = TRUE,
  public = list(
    # reg = NULL,
    cont.sVar.flag = logical(),
    bw.j = numeric(),
    glmfitclass = "glmS3", # default glm fit class
    is.fitted = FALSE,
    bindat = NULL, # object of class BinDat that is used in fitting / prediction, never saved (need to be initialized with $new())

    initialize = function(reg, ...) {
      assert_that(is.flag(reg$useglm))
      if (!reg$useglm) self$glmfitclass <- "speedglmS3"
      # self$reg <- reg
      self$bindat <- BinDat$new(reg = reg, ...) # postponed adding data in BinDat until self$fit() is called
      class(self$bindat) <- c(class(self$bindat), self$glmfitclass)
      # can also use: self$bindat <- BinDat$new(glm = self$glmfitclass, ...)
      # or: self$bindat <- BinDat$new(self, ...) (passing self might get confusing)
      if (gvars$verbose) {
        print("Init BinOutModel:"); print(self$show())
      }

      # Get the bin width (interval length) for the current bin name self$getoutvarnm (for discretized continuous sA only):
      self$cont.sVar.flag <- self$getoutvarnm %in% names(reg$intrvls.width)
      if (self$cont.sVar.flag) {
        intrvl.idx <- which(names(reg$intrvls.width) %in% self$getoutvarnm)
        if (length(intrvl.idx) > 1) stop("non-unique names for intrvls.width in RegressionClass")
        self$bw.j <- reg$intrvls.width[intrvl.idx]
      } else {
        self$bw.j <- 1L
      }

      invisible(self)
    },

    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked
      self$bindat$newdata(newdata = data, ...) # populate bindat with X_mat & Y_vals
      # self$bindat$newdata(newdata = data, getoutcome = TRUE, ...) # populate bindat with X_mat & Y_vals
      private$m.fit <- logisfit(datsum_obj = self$bindat) # private$m.fit <- data_obj$logisfit or private$m.fit <- data_obj$logisfit() 
      # alternative 2 is to apply data_obj method / method that fits the model
      self$is.fitted <- TRUE

      # **********************************************************************
      # to save RAM space when doing many stacked regressions:
      # **********************************************************************
      # if (self$reg$pool_cont && length(self$reg$outvars_to_pool) > 1) {
      #   private$probAeqa <- self$bindat$logispredict.long(m.fit = private$m.fit)
      # } else {
      #   private$probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
      # }

      # **********************************************************************
      # NOTE THIS WILL BREAK SUBSTITUTION ESTIMATOR FOR Q.n and Q.n^*, but is necessary to save space for many stacked regressions:
      # **********************************************************************
      # self$bindat$emptydata  # Xmat in bindat is no longer needed, only subset, outvar & probA1 (probOut1)
      self$wipe.alldat
      # self$bindat$emptyY
      # self$bindat$emptySubset_idx
      # self$bindat$emptyN

      invisible(self)
    },

    # take fitted BinOutModel class object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("BinOutModel" %in% class(bin.out.model))
      private$m.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # rename to:
    # predictP_1 = function(newdata, ...)
    predict = function(newdata, ...) { # P(A^s[i]=1|W^s=w^s): uses private$m.fit to generate predictions for newdata
      assert_that(self$is.fitted) # vs. stopifnot(self$is.fitted)
      if (missing(newdata)) {
        stop("must provide newdata for BinOutModel$predict()")
        # ... Do nothing, predictions for fitted data are already saved
        # return(invisible(self))
      }
      # re-populate bindat with new X_mat:
      self$bindat$newdata(newdata = newdata, ...) 
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        private$probAeqa <- self$bindat$logispredict.long(m.fit = private$m.fit)
      } else {
        private$probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
      }
      self$bindat$emptydata  # Xmat in bindat is no longer needed, but subset, outvar & probA1 may be needed for private$probA1
      invisible(self)
    },

    # take BinOutModel class object that contains the predictions for P(A=1|sW) and save these predictions to itself
    copy.predict = function(bin.out.model) {
      assert_that("BinOutModel" %in% class(bin.out.model))
      assert_that(self$is.fitted)
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        private$probAeqa <- bin.out.model$getprobAeqa
      } else {
        private$probA1 <- bin.out.model$getprobA1
      }
    },


    # **************************************
    # TO DO; Redo predictAeqa so that it has no side-effects, i.e., eliminate field private$probAeqa and only returns the result, without saving it.
    # **************************************
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    # Invisibly returns cumm. prob P(sA=sa|sW=sw)
    # rename to:
    # predict.like.P_a = function(newdata)
    predictAeqa = function(obs.DatNet.sWsA, bw.j.sA_diff) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a's)
      assert_that(self$is.fitted)
      assert_that(!missing(obs.DatNet.sWsA))
      self$bindat$newdata(newdata = obs.DatNet.sWsA) # populate bindat with new design matrix covars X_mat
      assert_that(is.logical(self$getsubset))
      n <- obs.DatNet.sWsA$nobs

      # obtain predictions (likelihood) for response on fitted data (from long pooled regression):
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        probAeqa <- self$bindat$logispredict.long(m.fit = private$m.fit) # overwrite probA1 with new predictions:

      } else {
        # get predictions for P(sA[j]=1|sW=newdata) from newdata:
        probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
        indA <- obs.DatNet.sWsA$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
        assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers

        probAeqa <- rep.int(1L, n) # for missing, the likelihood is always set to P(A = a) = 1.
        assert_that(!any(is.na(probA1[self$getsubset]))) # check that predictions P(A=1 | dmat) exist for all obs.
        probA1 <- probA1[self$getsubset]
        # discrete version for the joint density:
        probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
        # continuous version for the joint density:
        # probAeqa[self$getsubset] <- (probA1^indA) * exp(-probA1)^(1 - indA)

        # Alternative intergrating the last hazard chunk up to x:
        # difference of sA value and its left most bin cutoff: x - b_{j-1}
        if (!missing(bw.j.sA_diff)) {
          # + integrating the constant hazard all the way up to value of each sa:
          # probAeqa[self$getsubset] <- probAeqa[self$getsubset] * (1 - bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
          # cont. version of above:
          probAeqa[self$getsubset] <- probAeqa[self$getsubset] * exp(-bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
        }
      }
      # private$probAeqa <- probAeqa # disabling internal saving of probAeqa
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      # self$bindat$emptydata  # wipe out prediction data after getting the likelihood
      # self$bindat$emptyY
      # self$bindat$emptySubset_idx
      # private$probA1 <- NULL
      return(probAeqa)
    },

    show = function() {self$bindat$show()}

  ),
  active = list(
    wipe.alldat = function() {
      private$probA1 <- NULL
      private$probAeqa <- NULL
      self$bindat$emptydata
      self$bindat$emptyY
      self$bindat$emptySubset_idx
      self$bindat$emptyN
      return(self)
    },
    getfit = function() { private$m.fit },
    getprobA1 = function() { private$probA1 },
    getprobAeqa = function() { private$probAeqa },
    getsubset = function() { self$bindat$subset_idx },
    getoutvarval = function() { self$bindat$getY },
    getoutvarnm = function() { self$bindat$outvar }
  ),
  private = list(
    m.fit = list(),   # the model fit (coefficients)
    probA1 = NULL,    # Predicted probA^s=1 conditional on X_mat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on X_mat
  )
)