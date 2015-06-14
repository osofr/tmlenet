library(assertthat)
library(speedglm)
library(R6)
`%+%` <- function(a, b) paste0(a, b)  # cat function

NewSummaryModel = function(reg, ...) { UseMethod("NewSummaryModel") } # Generic S3 constructor for the summary model classes
NewSummaryModel.binary = function(reg, ...) {  # Summary model constructor for binary outcome sA[j]
  print("BinarySummaryModel constructor called...")
  BinarySummaryModel$new(glm = FALSE, reg = reg, ...) # fit a model with new object BinarySummaryModel
  # BinarySummaryModel$new(glm = FALSE, reg = reg, data = data, subset = subset) # fit a model with new object BinarySummaryModel  
}
NewSummaryModel.cat = function(reg, ...) { # Summary model constructor for categorical outcome sA[j]
  print("CatSummaryModel constructor called...")
  CatSummaryModel$new(...)
}
NewSummaryModel.contin = function(reg, ...) { # Summary model constructor for continuous outcome sA[j]
  print("ContinSummaryModel constructor called...")
  nbins <- 10L  # no. of bins for contin sA[j]
  # DEFINE THE NUMBER OF REGRESSIONS TO RUN, 
  # NEW OUTVAR var names, etc
  new.sA_class <- rep_len("binary", nbins)
  new.sA_nms <- "BIN_"%+%(1:nbins)%+%"_"%+%reg$outvar
  new.sW_nms <- reg$predvars

  new.subsets_chr <- lapply(new.sA_nms, function(var) {"!misfun("%+%var%+%")"})  # subsetting by !gvars$misval on sA:
  new.subset <- lapply(new.subsets_chr, function(subset_chr) {      
                                            subset_expr <- try(parse(text=subset_chr)[[1]]) # parses chr into a call
                                            if(inherits(subset_expr, "try-error")) stop("can't parse the subset formula", call.=FALSE)
                                            subset_expr
                                          })

  ContinSummaryModel$new(sA_class = new.sA_class, sA_nms = new.sA_nms, sW_nms = new.sW_nms, subset = new.subset, ...)
}

NewSummaryModel.continPOOL = function(reg, ...) { # Summary model constructor for continuous outcome sA[j]
  print("ContinSummaryModelPool constructor called...")
  # DEFINE THE NUMBER OF REGRESSIONS TO RUN, 
  # NEW OUTVAR var names, etc
  reg$outvar
  nbins <- 10L  # no. of bins for contin sA[j]
  new.sA_class <- rep_len("binary", nbins)
  new.sA_nms <- "BIN_"%+%reg$outvar
  new.sW_nms <- reg$predvars
  new.subset <- reg$subset
  # (sA_class, sA_nms, sW_nms = new.sW_nms, subset, ...)
  ContinSummaryModel$new(...)
}


## ---------------------------------------------------------------------
# Class for defining, managing, fitting and returning the likelihood P(sA = sa | sW = sw) under g_star or g_0;
# Accepts (1) data (data.frame) for (sA,sW), (2) newdata (data.frame) for prediction, (3) obsdat.sA (matrix) for sa values;
# Defines and manages the factorization of the joint P(sA = sa | ... ) into reg models sA[j] ~ \bar{sA[j-1]} + sW;
# Figures out reg mdel factorization based on name ordering in (sA_nms, sW_nms);
# Evaluates subset_exprs in the envirs of data and newdata data.frames
# Calls BinarySummaryModel$new, assumes each sA[j] is binary in reg (sA[j] ~ \bar{sA[j-1]} + sW);
## ---------------------------------------------------------------------
# TO DO 1: if (missing(sA_class)): init to a default (binary) vector of classes
# TO DO 2: if (missing(sA_class)): automatic class detection?
# TO DO 3: S3 dispatch for fit/predict based on sA_class[k_i] (reg$outvar.class)

#' @title Class for defining, holding and fitting collections of summary measure models P(sA[j]|sW,\bar{sA}[j])
#' @docType class
#' @format An R6 class object.
#' @name SummariesModel
#' @details Following fields are created during initialization
#' \itemize{
#' \item{n_regs} ...
#' \item{nodes} ...
#' \item{regs_list} ...
#' \item{sA_nms} ...
#' \item{sW_nms} ...
#' }
#' More details about the class...
##' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
##' @export
SummariesModel <- R6Class(classname = "SummariesModel",
  # inherit = Abstract_BinarySummaryModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    n_regs = integer(),     # total no. of reg. models (logistic regressions)
    # nodes = list(),         # list of node names in the data (Anode, Wnodes, nFnode, etc..); NOT USED YET
    regs_list = list(),     # List of regressions (summary measures) used for fitting. Each regression in regs_list is names in list(outvar =, predvar =)
    # subset_regs = list(),   # list of subsetting exprs (one logical expr per reg); NOT USED ANYMORE
    sA_nms = character(),   # sA names
    sW_nms = character(),   # sW names
    # Kmax = integer(),       # max n of Friends in the network; NOT USED YET

    initialize = function(sA_class, sA_nms, sW_nms, subset, ...) {
    # initialize = function(Kmax, nodes, sA_class, sA_nms, sW_nms, subset, ...) {
      # assert_that(is.count(Kmax))
      # self$Kmax <- Kmax
      # self$nodes <- nodes
      assert_that(is.character(sA_nms))
      assert_that(is.character(sW_nms))
      self$n_regs <- n_regs <- length(sA_nms)  # Number of sep. logistic regressions, determines the fitting algorithm this object will use (sep.logist by k, pooled, etc):
      # self$n_regs <- n_regs <- (Kmax + 1)
      self$sA_nms <- sA_nms
      self$sW_nms <- sW_nms

      print("#----------------------------------------------------------------------------------");
      print("New SummariesModel object:"); 
      print("No. of regressions: " %+% self$n_regs)
      print("sA_classes: "); print(sA_class)
      print("sA_nms: "); print(sA_nms)
      print("sW_nms: "); print(sW_nms)
      # print("subsets: "); print(subset)
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
        sA_i_nm <- sA_nms[k_i] # A variable we are predicting
        covars_nms <- c(sA_nms[-c(k_i:n_regs)], sW_nms) # dependent covars

        reg <- list(outvar.class = sA_class[k_i], outvar = sA_i_nm, predvars = covars_nms, subset = subset[[k_i]])
        class(reg) <- c(class(reg), reg$outvar.class)

        # Declare new object for summary model P(sA[j]|\bar{sA}[j-1], sW}):
        PsAsW.model <- NewSummaryModel(reg = reg) 
        private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model)) # save full fitted BinarySummaryModel:
        names(private$PsAsW.models)[k_i] <- "P(sA|sW)."%+%k_i

        self$regs_list <- append(self$regs_list, list(reg)) # self$subset_regs <- append(self$subset_regs, subset[[k_i]])  # removed (added to reg list)
      }
      names(self$regs_list)[1:n_regs] <- "reg."%+%1:n_regs

      # names(self$subset_regs)[1:n_regs] <- "determ.reg."%+%1:n_regs
      invisible(self)
    },
    length = function(){ base::length(self$regs_list) },
    getfitted.pbins = function() { private$fitted.pbins },  # get all summaries objs
    getcumprodAeqa = function() { private$cumprodAeqa },  # get a vector of cum prod of P(sA[j]=a[j]|sW)
    fit = function(data) {  
      n <- nrow(data)
      for (k_i in seq_along(self$regs_list)) { # loop over all regressions in regs_list
        # replaced below with:
        private$PsAsW.models[[k_i]]$fit(data = data)
        # pbin <- BinarySummaryModel$new(glm = FALSE, reg = reg, data = data, subset = subset) # fit a model with new object BinarySummaryModel
        # pbin$fit(reg = reg, data = data, subset = subset) # fit a model with new object BinarySummaryModel  
        
        print("#----------------------------------------------------------------------------------");
        print("SummariesModel$fit(data), private$PsAsW.models[[k_i]]:");
        print(private$PsAsW.models[[k_i]]$show())
        print("self$regs_list[[k_i]]: "); str(self$regs_list[[k_i]]);
        print("PsAsW.models$datbin$getY: "); print(head(private$PsAsW.models[[k_i]]$datbin$getY))
        print("PsAsW.models[[k_i]]$getsubset: "); print(head(private$PsAsW.models[[k_i]]$getsubset))
        print("#----------------------------------------------------------------------------------");
      }
      invisible(self)
    },
    # use BinarySummaryModel objects (stored from prev call to fit) to predict P(A=1|..) based on new cY_i matrix
    predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
      if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
        return(invisible(self))
      }
      n <- nrow(newdata)
      for (k_i in seq_along(self$regs_list)) { # loop over all regressions in regs_list
        private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        # reg <- self$regs_list[[k_i]] # subset_expr <- self$subset_regs[[k_i]]
        # private$fitted.pbins[[k_i]]$predict(newdata = newdata, subset = subset)
      }
      invisible(self)
    },
    # use BinarySummaryModel objects (stored from prev call to fit()) to predict P(sA=obsdat.sA|sW)
    predictAeqa = function(obsdat.sA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      assert_that(is.matrix(obsdat.sA)) # check A: obsdat.sA is a matrix
      n <- nrow(obsdat.sA)
      cumprodAeqa <- rep_len(1, n)
      for (k_i in seq_along(self$regs_list)) { # loop over all regressions in regs_list
        cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(indA_i = obsdat.sA[, k_i])$getprobAeqa()
      }
      private$cumprodAeqa <- cumprodAeqa
      invisible(self)
    }
  ),

  active = list(
    actplhold = function() {}
    # x2 = function(value) {
    #   if (missing(value)) return(self$x * 2)
    #   else self$x <- value/2
    # },
    # rand = function() {rnorm(1)}
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
CatSummaryModel <- R6Class(classname = "CatSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE
)

ContinSummaryModel <- R6Class(classname = "ContinSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    # remains unchanged: initialize = function(Kmax, nodes, sA_nms, sW_nms, subset, ...) {      


# # NORMALIZE:
#       x <- (x - min(x)) / (max(x) - min(x))
# # 2A) turn x into an ordinal (0, 1, 2, 3, ...) based on interval defn 1: bin x by equal length intervals of 0-1:
#       intlen <- .1; # generates 10 intervals (intlen <- .2; intlen <- .5;)
#       intvec <- seq(from = 0, to = 1, by=intlen)
#       system.time(
#         i <- data.frame(x = x, xcat = findInterval(x = x, vec = intvec, rightmost.closed = TRUE))
#       )
#       # for 100K, intlen=.1:
#       # user  system elapsed 
#       # 0.016   0.001   0.026 
#       # length(unique(i$xcat)); hist(i$xcat)
#       # 2B) turn x into an ordinal (0, 1, 2, 3, ...) based on interval defn 2: bin x by mass (quantiles on 0-1 as probs)
#       system.time( {
#         quantvec <- quantile(x=x, probs = intvec)
#         i["xcatq"] <- findInterval(x = x, vec = quantvec, rightmost.closed = TRUE)    
#       })      
#       # length(unique(i$xcatq)); hist(i$xcatq)

# # 3) CREATE INDICATOR DUMMIES FOR ALL NON-BINARY SUMMARY MEASUREs in sA (do not touch sW)
# # Approach 1: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
      
#       # ***
#       # TO DO: Replace degenerate indicators with gvars$misval
#       # The output needs to be a df because of subset evaluation
#       # Is it possible to unify NetInd_i way of creating matrix (NA for missing friends) with this?
#       # ***

#       cats <- sort(unique(i$xcat))
#       bvarnmsnew <- "B"%+%"_"%+%cats[-length(cats)]
#       system.time({
#         dummies2 <- matrix(nrow=nrow(i), ncol=length(cats)-1)
#         for(level in cats[-length(cats)]) {
#             # dummies2[,level] <- as.integer(i$xcat <= level)
#           # rewrite of above (***NOT TESTED***):
#           subset_Bj <- i$xcat <= level
#           dummies2[subset_Bj, level] <- as.integer(subset_Bj)
#           dummies2[!subset_Bj, level] <- gvars$misval
#         }
#         colnames(dummies2) <- bvarnmsnew
#       })
#       head(dummies2)
#       class(dummies2)
#       class(dummies2[,1])
#       # for 100K numeric RVs with intlen <- .1;
#       #  user  system elapsed 
#       # 0.037   0.006   0.054
#       # 0.032   0.005   0.038 
# # Approach 2 (model.matrix): a bit slower than manual
#       # Creates B_j that jumps to 1 and then back to 0 (degenerate) and uses category 1 as ref (droped)
#       # First need to make a factor from xcat; adds a column of 1's (removed manually)
#       cats <- sort(unique(i$xcat))
#       bvarnmsnew <- "B"%+%"_"%+%cats[-1]
#       system.time({
#         i$xcatf <- factor(i$xcat);
#         dummies1 <- model.matrix(~i$xcatf);
#         saveattrs <- attributes(dummies1)
#         bvarnms <- colnames(dummies1)[-1]
#         colnames(dummies1)[-1] <- bvarnmsnew
#         # i <- data.frame(i, dummies1[,-1])
#       })
#       # Approach 2 manual. Creates B_j that jumps to 1 and then back to 0 (degenerate) and exclude the reference category (last)
#       cats <- sort(unique(i$xcat))
#       bvarnmsnew <- "B"%+%"_"%+%cats[-length(cats)]
#       system.time({
#         dummies2 <- matrix(nrow=nrow(i), ncol=length(cats)-1)
#         for(level in cats[-length(cats)]) {
#           dummies2[,level] <- as.integer(i$xcat == level)
#         }
#         colnames(dummies2) <- bvarnmsnew
#       })




    # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      n <- nrow(data)
      # ...
      # transform data (either define binary columns (BinsA[1],...,BinsA[M]) or pool the binaries into one colum dataset)
      # ...
      data_t <- data
      super$fit(data_t) # call parent fit method
      invisible(self)
    },
    # use BinarySummaryModel objects (stored from prev call to fit) to predict P(A=1|..) based on new cY_i matrix
    predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
      if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
        return(invisible(self))
      }
      # ...
      # transform newdata (either define binary columns (BinsA[1],...,BinsA[M]) or pool the binaries into one colum dataset)
      # ...
      newdata_t <- newdata
      super$predict(newdata_t)
      invisible(self)
    },
    # use BinarySummaryModel objects (stored from prev call to fit()) to predict P(sA=obsdat.sA|sW)
    predictAeqa = function(obsdat.sA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      assert_that(is.matrix(obsdat.sA)) # check A: obsdat.sA is a matrix
      n <- nrow(obsdat.sA)
      # ...
      # transform obsdat.sA (either define binary columns (BinsA[1],...,BinsA[M]) or pool the binaries into one colum dataset)
      # ...
      obsdat.sA_t <- obsdat.sA
      obsdat.sA_t
      invisible(self)
    }
  ),

  active = list(
    actplhold = function() {}
    # x2 = function(value) {
    #   if (missing(value)) return(self$x * 2)
    #   else self$x <- value/2
    # },
    # rand = function() {rnorm(1)}
  ),

  private = list(
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)
