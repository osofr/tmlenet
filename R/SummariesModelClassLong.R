#----------------------------------------------------------------------------------
# Classes that control modelling of the sequential G-computation for multiple point-treatments
#----------------------------------------------------------------------------------
#' @importFrom assertthat assert_that
is.SGCompRegClass <- function(obj) "SGCompRegClass"%in%class(obj)

## ---------------------------------------------------------------------
#' @export
# FOR NOW ASSUMING STRUCTURE O=(W,A(0),C(0),L(0),A(1),C(1),L(1),...Y)
SGCompRegClass <- R6Class("SGCompRegClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    outvar = character(),          # the outcome variable name (Ynode)
    outvar.class = character(),

    Qforms = character(),
    n_regs = integer(),
    # n_timepts = integer(),       # number of time points
    nodes = list(),
    predvars = list(),             # list of predictors for each regression from Qforms

    # not used for now:
    Anodes = list(),               # list of treatment var names, by timepoint
    Cnodes = list(),               # list of censoring var names, by timepoint
    Lnodes = list(),               # list of time-varying confounder names, by timepoint
    Wnodes = character(),          # character vector of baseline covariate names
    subset = NULL,                 # subset expression (later evaluated to logical vector in the envir of the data)

    ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    pool_cont = FALSE,             # NOT USED IN G-COMP
    outvars_to_pool = NULL,        # NOT USED IN G-COMP

    useglm = logical(),            # TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
    parfit = logical(),            # TRUE for fitting binary regressions in parallel
    pool = logical(),

    # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
    # form = NULL,                 # (NOT IMPLEMENTED) reg formula, if provided run using the usual glm / speedglm functions
    initialize = function(outvar,
                          nodes,
                          Qforms,
                          # n_timepts,
                          # Anodes, Cnodes, Lnodes, Wnodes,
                          # predvars,
                          # subset,
                          # ReplMisVal0 = TRUE, # Needed to add ReplMisVal0 = TRUE for case sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise:
                          useglm = getopt("useglm"),
                          parfit = getopt("parfit"),
                          pool = FALSE,
                          family = "quasibinomial") {


      self$outvar <- outvar
      # always doing univariate regression on the outcome (even if continuous):
      self$outvar.class <- gvars$sVartypes$bin
      self$nodes <- nodes
      self$Qforms <- Qforms
      self$predvars <- lapply(Qforms, RhsVars)
      # self$ReplMisVal0 <- ReplMisVal0
      self$useglm <- useglm
      self$parfit <- parfit
      self$pool <- pool
      self$n_regs <- n_regs <- length(Qforms)

      # if (!missing(subset)) {
      #   self$subset <- subset
      #   if (length(subset) < n_regs) {
      #     self$subset <- rep_len(subset, n_regs)
      #   } else if (length(subset) > n_regs) {
      #     # ... TO FINISH ...
      #     if (!is.logical(subset)) stop("not implemented")
      #     # increase n_regs to all combinations of (n_regs x subset)
      #   }
      # } else {
      #   self$subset <- rep_len(list(TRUE), n_regs)
      # }
    },

    # take the clone of a parent RegressionClass (reg) for length(self$outvar) regressions
    # and set self to a single univariate k_i regression for outcome self$outvar[[k_i]]
    ChangeManyToOneRegresssion = function(k_i, reg) {
      assert_that(!missing(k_i))
      if (missing(reg)) stop("reg must be specified along with k_i")
      assert_that(is.count(k_i))
      assert_that(k_i <= reg$n_regs)
      self$n_regs <- 1
      self$outvar.class <- gvars$sVartypes$bin
      # self$outvar.class <- reg$outvar.class[[k_i]] # Class of the outcome var: binary, categorical, continuous:

      # might change to time-point specific intermediate outcome names:
      self$outvar <- "Q.kplus1"
      # self$outvar <- reg$outvar[[k_i]] # An outcome variable that is being modeled:

      self$predvars <- reg$predvars[[k_i]] # Regression covars (predictors):

      # the subset is a list when RegressionClass specifies several regression models at once,
      # obtain the appropriate subset for this regression k_i and set it to self
      if (is.list(reg$subset)) {
        self$subset <- reg$subset[[k_i]]
      }
      self$S3class <- self$outvar.class # Set class on self for S3 dispatch...
      return(invisible(self))
    },

    resetS3class = function() class(self) <- c("RegressionClass", "R6")

  ),

  active = list(
    # For S3 dispatch on newsummarymodel():
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
      list(outvar = self$outvar,
          predvars = self$predvars,
          subset = self$subset)
    }
  )
)

#' @export
SGcompSummariesModel <- R6Class(classname = "SGcompSummariesModel",
  portable = TRUE,
  class = TRUE,
  public = list(
    # reg = NULL,
    outvar = character(),      # the name of the continous outcome var (sA[j])
    n_regs = integer(),        # total no. of reg. models (logistic regressions)
    parfit_allowed = FALSE,    # allow parallel fit of multivar outvar when 1) reg$parfit = TRUE & 2) all.outvar.bin = TRUE
    initialize = function(reg, ...) {
      assert_that(is.SGCompRegClass(reg))
      self$n_regs <- reg$n_regs # Number of sep. sequential regressions to run
      self$outvar <- reg$outvar

      # if (gvars$verbose) {
        print("#----------------------------------------------------------------------------------")
        print("New instance of SGcompSummariesModel:")
        print("#----------------------------------------------------------------------------------")
        print("Outcomes: " %+% paste(reg$outvar, collapse = ", "))
        # print("nodes: " %+% paste(unlist(reg$nodes), collapse = ", "))
        # print("predvars: " %+% paste(unlist(reg$predvars), collapse = ", "))
        print("No. of regressions: " %+% self$n_regs)
        # print("Anodes: " %+% paste(unlist(reg$Anodes), collapse = ", "))
        # print("Cnodes: " %+% paste(unlist(reg$Cnodes), collapse = ", "))
        # print("Lnodes: " %+% paste(unlist(reg$Lnodes), collapse = ", "))
        # print("Wnodes: " %+% paste(reg$Wnodes, collapse = ", "))
        print("#----------------------------------------------------------------------------------")
      # }

      # factorize E[Y_d|L(0)] into GCOMP sequential regressions, by number of time points:
      # start from the very end (regressing Y on all past, all the way up to first regression on W)
      for (k_i in self$n_regs:1) {
        reg_i <- reg$clone()
        reg_i$ChangeManyToOneRegresssion(k_i, reg)
        # Calling the constructor for the summary model P(sA[j]|\bar{sA}[j-1], sW}), dispatching on reg_i class
        PsAsW.model <- newsummarymodel(reg = reg_i, ...)
        private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
        names(private$PsAsW.models)[length(private$PsAsW.models)] <- "P(sA|sW)."%+%k_i
      }
      invisible(self)
    },

    # Sequential G-COMP regression
    # TO DO: NEED TO ADD newdata ARG TO THIS WHICH WILL CONTAIN COUNTERFACTUAL TREATMENT ASSIGNMENTS
    fit = function(data, newdata) {
    # fit = function(data, newdata) {
      assert_that(is.DatNet.sWsA(data))
      # if (gvars$verbose) {
        print("fitting sequential G-comp: " %+% self$outvar)
      # }

      # Initialize Q.kplus1 by adding outcome to data$Q.kplus1:
      data$initQ.kplus1(n_regs = self$n_regs)

      # sequential regression loop for Gcomp, alternate fit & predict:
      for (k_i in seq_along(private$PsAsW.models)) {

        # 1. run one iteration of the seq. G-comp regression:
        # private$PsAsW.models[[k_i]]$fit(data = data, getoutvar = ifelse(k_i==1, TRUE, "Q.kplus1"))
        private$PsAsW.models[[k_i]]$fit(data = data, getoutvar = "Q.kplus1")

        # 2. set Anodes to specific regimens (gstar) so that we can predict the intermediate counterfactual E[Y_d|...A(t)=a(t),L(t)):
        # ... NOT IMPLEMENTED ...

        # 3. predict intermediate Y_d, based on Anodes set to gstar, this predictions will be used as the outcome during next regression:
        Q.kplus1 <- private$PsAsW.models[[k_i]]$predict(newdata = newdata)$getprobA1
        # Q.kplus1 <- private$PsAsW.models[[k_i]]$predict(newdata = data)$getprobA1

        # 4. record Q.kplus1, so that its available as an outcome to next regression:
        data$addQ.kplus1(name = "Q.kplus1", iter = k_i, val = Q.kplus1)
        print("Q.kplus1 fit for sequential regression: " %+% k_i); print(head(Q.kplus1))
      }

      # clean up:
      # data$emptydat.bin.sVar # wiping out binirized mat in data DatNet.sWsA object...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...

      invisible(self)
    },

    # ******* NOT USED *******
    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) {
        stop("must provide newdata")
      }
      assert_that(is.DatNet.sWsA(newdata))
      # serial loop over all regressions in PsAsW.models:
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$predict(newdata = newdata)
      }
      invisible(self)
    },

    # getcumprodAeqa = function() { private$cumprodAeqa },  # get joint prob as a vector of the cumulative prod over j for P(sA[j]=a[j]|sW)
    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Uses daughter objects (stored from prev call to fit()) to get predictions for P(sA=obsdat.sA|sW=sw)
    # Invisibly returns the joint probability P(sA=sa|sW=sw), also saves it as a private field "cumprodAeqa"
    # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    # predictAeqa = function(newdata, ...) {
    #   assert_that(!missing(newdata))
    #   assert_that(is.DatNet.sWsA(newdata))
    #   n <- newdata$nobs
    #   cumprodAeqa <- rep.int(1, n)
    #   # loop over all regressions in PsAsW.models:
    #   for (k_i in seq_along(private$PsAsW.models)) {
    #       cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, ...)
    #   }
    #   private$cumprodAeqa <- cumprodAeqa
    #   return(cumprodAeqa)
    # },

    length = function(){ base::length(private$PsAsW.models) },
    getPsAsW.models = function() { private$PsAsW.models }  # get all summary model objects (one model object per outcome var sA[j])
  ),

  active = list(
    # recursively call all saved daughter model fits and wipe out any traces of saved data:
    wipe.alldat = function() {
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$wipe.alldat
      }
      return(self)
    }
  ),

  private = list(
    deep_clone = function(name, value) {
      # if value is is an environment, quick way to copy:
      # list2env(as.list.environment(value, all.names = TRUE), parent = emptyenv())
      # if a list of R6 objects, make a deep copy of each:
      if (name == "PsAsW.models") {
        lapply(value, function(PsAsW.model) PsAsW.model$clone(deep=TRUE))
      # to check the value is an R6 object:
      } else if (inherits(value, "R6")) {
        value$clone(deep=TRUE)
      } else {
        value # For all other fields, just return the value
      }
    },
    PsAsW.models = list(),
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)


# -------------------------------------------------------------------------------------------
# same code in ContinSummaryModel$new and CategorSummaryModel$new replaced with outside function:
# Define subset evaluation for new bins:
# ******************************************************
# NOTE: Subsetting by var name only (which automatically evaluates as !gvars$misval(var)) for speed & memory efficiency
# ******************************************************
# -------------------------------------------------------------------------------------------
# def_regs_subset <- function(self) {
#   bin_regs <- self$reg$clone() # instead of defining new RegressionClass now cloning parent reg object and then ADDING new SETTINGS
#   bin_regs$reg_hazard <- TRUE # don't add degenerate bins as predictors in each binary regression
#   if (!self$reg$pool_cont) {
#     add.oldsubset <- TRUE
#     new.subsets <- lapply(self$reg$bin_nms,
#                               function(var) {
#                                 res <- var
#                                 if (add.oldsubset) res <- c(res, self$reg$subset)
#                                 res
#                               })

#     new.sAclass <- as.list(rep_len(gvars$sVartypes$bin, self$reg$nbins))
#     names(new.sAclass) <- self$reg$bin_nms
#     bin_regs$ChangeOneToManyRegresssions(regs_list = list(outvar.class = new.sAclass,
#                                                           outvar = self$reg$bin_nms,
#                                                           predvars = self$reg$predvars,
#                                                           subset = new.subsets))
#   # Same but when pooling across bin indicators:
#   } else {
#     bin_regs$outvar.class <- gvars$sVartypes$bin
#     bin_regs$outvar <- self$outvar
#     bin_regs$outvars_to_pool <- self$reg$bin_nms
#     if (gvars$verbose)  {
#       print("pooled bin_regs$outvar: "); print(bin_regs$outvar)
#       print("bin_regs$outvars_to_pool: "); print(bin_regs$outvars_to_pool)
#       print("bin_regs$subset: "); print(bin_regs$subset)
#     }
#   }
#   bin_regs$resetS3class()
#   return(bin_regs)
# }

# # -------------------------------------------------------------------------------------------
# #' @export
# ContinSummaryModel <- R6Class(classname = "ContinSummaryModel",
#   inherit = SummariesModel,
#   portable = TRUE,
#   class = TRUE,
#   public = list(
#     reg = NULL,
#     outvar = character(),     # the name of the continous outcome var (sA[j])
#     intrvls = NULL,
#     intrvls.width = NULL,
#     bin_weights = NULL,
#     # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
#     initialize = function(reg, DatNet.sWsA.g0, DatNet.sWsA.gstar, ...) {
#       self$reg <- reg
#       self$outvar <- reg$outvar
#       self$reg$bin_nms <- DatNet.sWsA.g0$bin.nms.sVar(reg$outvar, self$reg$nbins)
#       # Save bin widths in reg class (naming the vector entries by bin names):
#       if (gvars$verbose)  {
#         print("ContinSummaryModel outcome: "%+%self$outvar)
#       }
#       bin_regs <- def_regs_subset(self = self)
#       super$initialize(reg = bin_regs, ...)
#     },

#     # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
#     # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
#     fit = function(data) {
#       assert_that(is.DatNet.sWsA(data))
#       # Binirizes & saves binned matrix inside DatNet.sWsA
#       data$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
#       if (gvars$verbose) {
#         print("performing fitting for continuous outcome: " %+% self$outvar)
#         print("freq counts by bin for continuous outcome: "); print(table(data$ord.sVar))
#         print("binned dataset: "); print(head(cbind(data$ord.sVar, data$dat.bin.sVar), 5))
#       }
#       super$fit(data) # call the parent class fit method
#       if (gvars$verbose) message("fit for outcome " %+% self$outvar %+% " succeeded...")
#       data$emptydat.bin.sVar # wiping out binirized mat in data DatNet.sWsA object...
#       self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
#       invisible(self)
#     },

#     # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
#     predict = function(newdata) {
#       if (missing(newdata)) {
#         stop("must provide newdata")
#       }
#       assert_that(is.DatNet.sWsA(newdata))
#       if (gvars$verbose) print("performing prediction for continuous outcome: " %+% self$outvar)
#       # mat_bin doesn't need to be saved (even though its invisibly returned); mat_bin is automatically saved in datnet.sW.sA - a potentially dangerous side-effect!!!
#       newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
#       super$predict(newdata)
#       newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DatNet.sWsA object...
#       invisible(self)
#     },

#     # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
#     # Invisibly return cumm. prob P(sA=sa|sW=sw)
#     predictAeqa = function(newdata) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
#       assert_that(is.DatNet.sWsA(newdata))
#       newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
#       if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
#       bws <- newdata$get.sVar.bw(name.sVar = self$outvar, intervals = self$intrvls)
#       self$bin_weights <- (1 / bws) # weight based on 1 / (sVar bin widths)
#       # Option 1: ADJUST FINAL PROB by bw.j TO OBTAIN density at a point f(sa|sw) = P(sA=sa|sW=sw):
#       cumprodAeqa <- super$predictAeqa(newdata = newdata) * self$bin_weights
#       newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
#       self$bin_weights <- NULL # wiping out self$bin_weights...
#       self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
#       private$cumprodAeqa <- cumprodAeqa
#       return(cumprodAeqa)
#     }
#   ),
#   active = list(
#     cats = function() {seq_len(self$reg$nbins)}
#   )
# )

