

# **********
# TO DO (Continuous summary measures)
# x) See how to generalize to pooled fits, k-specific fits, etc (use subset definitions + ?)
# x) Add function to convert data_mtx to long format for SummaryM.pool
# x) For SummaryM.pool how to define predict function?

# **********
# TO DO (ContinSummaryModel): 
# * MOVE binarize() fun to DatNet and make it part of DatNet experience.
# * Need to handle gvars$misval values for self$outvar (contin) when transforming to bins (in fit, predict, predictAeqa)...
# * For x_cat the process is identical, except that normalize, define.intervals & make.ordinal is skipped.
# * NEED TO TEST that make.bins_mtx_1 will do the right thing when x_cat is (0, ..., ncats) instead of (1, ..., ncats) (IT SHOULD WORK)

# **********
# TO DO (NewSummaryModel.contin):
# * Need to find a way to pass nbins/bin_bymass/add.oldsubset for each contin/cat sA[j] from DatNet class...
# * These args need to be generated in some automated and consistent way for all SummariesModel$new()
# * Remove subset_expr definition from here, where is a more approapriate location?


RegressionClass <- R6Class("RegressionClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    pool.cont = FALSE,          # pool observations (long format) for binarized continous outcomes from different bins
    useglm = TRUE,              # TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
    outvar.class = character(), # class of the outcome var: bin / cont / cat
    outvar = character(),       # regression outcome variable name
    predvars = character(),     # vector of regression covariate names (predictors)
    subset = NULL,              # subset expression (evaluated to logical subset vector in the environment of the data)
    initialize = function(outvar.class, outvar, predvars, subset) {
      self$outvar.class <- outvar.class
      class(self) <- c(self$outvar.class, class(self))
      self$outvar <- outvar
      self$predvars <- predvars
      self$subset <- subset
    }
  ),
  active = list(
    reg = function() { list(outvar.class = self$outvar.class,
                            outvar = self$outvar,
                            predvars = self$predvars,
                            subset = self$subset) }
  )
)

NewSummaryModel = function(reg, datnet.sA, ...) { UseMethod("NewSummaryModel") } # Generic S3 constructor for the summary model classes

NewSummaryModel.contin = function(reg, datnet.sA, ...) { # Summary model constructor for continuous outcome sA[j]
  print("ContinSummaryModel constructor called...")
  ContinSummaryModel$new(reg = reg, datnet.sA = datnet.sA, ...)
  # Alternative name: ContOutModel
}

# NewSummaryModel.cat = function(reg, ...) { # Summary model constructor for categorical outcome sA[j]
#   print("CatSummaryModel constructor called...")
#   CatSummaryModel$new(...)
# }

NewSummaryModel.binary = function(reg, ...) {  # Summary model constructor for binary outcome sA[j]
  print("BinarySummaryModel constructor called...")
  # BinarySummaryModel$new(glm = TRUE, reg = reg, ...) # fit a model with new object BinarySummaryModel
  BinarySummaryModel$new(glm = FALSE, reg = reg, ...) # fit a model with new object BinarySummaryModel
  # Alternative name: BinOutModel
}


## ---------------------------------------------------------------------
# Class for defining, managing, fitting and returning the likelihood P(sA = sa | sW = sw) under g_star or g_0;
# Accepts (1) data (data.frame) for (sA,sW), (2) newdata (data.frame) for prediction, (3) obsdat.sA (matrix) for sa values;
# Defines and manages the factorization of the joint P(sA = sa | ... ) into reg models sA[j] ~ \bar{sA[j-1]} + sW;
# Figures out reg mdel factorization based on name ordering in (sA_nms, sW_nms);
# Evaluates subset_exprs in the envirs of data and newdata data.frames
# Calls BinarySummaryModel$new, assumes each sA[j] is binary in reg (sA[j] ~ \bar{sA[j-1]} + sW);
## ---------------------------------------------------------------------

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
    regs_list = list(),     # List of regressions (summary measures) used for fitting. Each regression in regs_list is names in list(outvar =, predvar =)
    sA_nms = character(),   # sA names
    sW_nms = character(),   # sW names

    initialize = function(sA_class, sA_nms, sW_nms, subset, ...) {
      # assert_that(is.character(sA_nms))
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
        self$regs_list <- append(self$regs_list, list(reg))
      }
      names(self$regs_list)[1:n_regs] <- "reg."%+%1:n_regs
      invisible(self)
    },
    show = function() {str(self$regs_list)},
    length = function(){ base::length(self$regs_list) },
    getfitted.pbins = function() { private$fitted.pbins },  # get all summaries objs
    getcumprodAeqa = function() { private$cumprodAeqa },  # get a vector of cum prod of P(sA[j]=a[j]|sW)
    fit = function(data) {  
      # n <- nrow(data)
      for (k_i in seq_along(self$regs_list)) { # loop over all regressions in regs_list        
        private$PsAsW.models[[k_i]]$fit(data = data) # below is replaced with this
        # print("#----------------------------------------------------------------------------------");
        # print("SummariesModel$fit(data), private$PsAsW.models[[k_i]]:");
        # print("k_i: "%+%k_i);
        # print("self$regs_list[[k_i]]: "); str(self$regs_list[[k_i]]);
        # print(private$PsAsW.models[[k_i]]$show())
        # print("PsAsW.models$datbin$getY: "); print(head(private$PsAsW.models[[k_i]]$datbin$getY))
        # print("PsAsW.models[[k_i]]$getsubset: "); print(head(private$PsAsW.models[[k_i]]$getsubset))
        # print("#----------------------------------------------------------------------------------");
      }
      invisible(self)
    },
    # use BinarySummaryModel objects (stored from prev call to fit) to predict P(A=1|..) based on new cY_i matrix
    predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
      if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
        return(invisible(self))
      }
      for (k_i in seq_along(self$regs_list)) { # loop over all regressions in regs_list
        private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        # private$fitted.pbins[[k_i]]$predict(newdata = newdata, subset = subset)
      }
      invisible(self)
    },
    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Use daughter objects (stored from prev call to fit()) to run predict on P(sA=obsdat.sA|sW)
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(obsdat.sA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      assert_that(is.matrix(obsdat.sA)) # check A: obsdat.sA is a matrix
      n <- nrow(obsdat.sA)
      cumprodAeqa <- rep_len(1, n)
      for (k_i in seq_along(self$regs_list)) { # loop over all regressions in regs_list
        cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(obsdat.sA = obsdat.sA[, k_i])
        print("getting likelihood for: "%+%self$regs_list[[k_i]]$outvar)
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
    datnet.sA = NULL,
    outvar = character(),     # the name of the continous outcome var (sA[j])
    nbins = integer(),        # no. of bins for discretizing contin/cat sA[j]
    # bin_bymass = TRUE,      # by default the bin cutoffs are determined using the quantiles of normalized sA[j]
    bin_intrvls = NULL,       # intervals (cut-off points) that define bins, calculated at $fit(data) stage
    bin_nms = NULL,

    initialize = function(reg, datnet.sA, ...) {  # define settings for fitting contin sA and then call $new for super class (SummariesModel)         
      # Define the number of bins (no. of regs to run), new outvar var names, prevars remain unchanged:
      self$datnet.sA <- datnet.sA
      self$outvar <- reg$outvar
      self$nbins <- datnet.sA$nbins.sVar(reg$outvar)
      self$bin_intrvls <- datnet.sA$get.sVar.intrvls(reg$outvar)

      self$bin_nms <- new.sA_nms <- datnet.sA$bin.nms.sVar(reg$outvar) #old: new.sA_nms <- reg$outvar%+%"_"%+%"B."%+%(1:cont.reg_par$nbins)
      new.sA_class <- as.list(rep_len(gvars$sVartypes$bin, self$nbins)); names(new.sA_class) <- new.sA_nms
      bin.m.params <- list(sA_class = new.sA_class, sA_nms = new.sA_nms, sW_nms = reg$predvars)

      # *) new.subsets:
      # TO DO: Put this in a separate function (with var as arg + additional args) +
      # move new.subsets def into another location (inside ContinSummaryModel$new()?)
      add.oldsubset <- TRUE
      new.subsets_chr <- lapply(new.sA_nms, function(var) {
                                              newsub <- "!misfun("%+%var%+%")"
                                              if (add.oldsubset) {
                                                newsub <- newsub %+% " & (" %+% deparse(reg$subset) %+% ")"
                                              }
                                              newsub
                                            })
      new.subsets <- lapply(new.subsets_chr, function(subset_chr) {      
                                                subset_expr <- try(parse(text=subset_chr)[[1]])
                                                if(inherits(subset_expr, "try-error")) stop("can't parse the subset formula", call.=FALSE)
                                                subset_expr
                                              })
      bin.m.params <- append(bin.m.params, list(subset = new.subsets))


      # Combine list of params binparams from datnet.sA with additinal params passed in ...:
      addl_params <- list(...)
      bin.m.params <- append(bin.m.params, addl_params)
      parnames <- names(bin.m.params)
      if (length(bin.m.params) != 0 && (is.null(parnames) || any(parnames==""))) {
        stop("please specify name for each attribute")
      }

      print("new in continuous summary...")
      print("contin sA name: "%+%self$outvar);
      print("contin sA nbins"); print(self$nbins)
      print("contin sA bin_intrvls"); print(self$bin_intrvls)
      print("contin binned sA names: "); print(new.sA_nms)
      print("Contin. new.subsets: "); str(new.subsets)
      # print("contin bin_bymass? "%+%self$bin_bymass)

      do.call(super$initialize, bin.m.params)
      # super$initialize(...) # call the parent class contructor
      # super$initialize(sA_class = new.sA_class, sA_nms = new.sA_nms, sW_nms = new.sW_nms, subset = new.subsets, ...) # call the parent class contructor
    },

    # binirize = function(contin.sAj) { # Transforms continous outcome sA[j] into matrix of discretized bin columns (sA[j] -> BinsA[1], ..., BinsA[M])
    #   contin.sAj <- normalize(x = contin.sAj) # norm 0-1, this step is optional
    #   if (is.null(self$bin_intrvls)) { # The intervals (cut-off points) defining bins need to be calculated only once and then saved:
    #     self$bin_intrvls <- define.intervals(x = contin.sAj, nbins = self$nbins, bin_bymass = self$bin_bymass) # define cut-off points
    #   }
    #   ord.sAj <- make.ordinal(x = contin.sAj, intervals = self$bin_intrvls) # transform data (either define binary columns (BinsA[1],...,BinsA[M]) or pool the binaries into one colum dataset?)
    #   # print("ord.sAj: "); print(ord.sAj)
    #   bin.sAj_mat <- make.bins_mtx_1(x.ordinal = ord.sAj, self = self) # print("bin.sAj_mat: "); print(head(bin.sAj_mat))
    #   make.bins_mtx_1(x.ordinal, nbins = self$nbins, bin.nms = ) { # Make dummy indicators for continuous x (sA[j])
    #   return(bin.sAj_mat)
    # },

    # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      print("fit in continuous summary...")
      # data is of class datnet.sA
      # data$datnetW$dat.sVar -> returns the matrix of sW
      # data$datnetA$dat.sVar -> returns a matrix of sA
      # data$get.df.sW.sA -> returns a data.frame of sWsA

      # TO DO: SORT OUT SAVING / SUBSETING / RETURNING mat_bin in data 
      # Evaluates subsets & return the correct data.frame to DatBin
      mat_bin <- data$binirize.sVar(self$outvar)
      # mat_bin <- data$datnetA$binirize.sVar(self$outvar)
      # print("mat_bin"); print(head(mat_bin))
      # print("data$bin_mat.sVar"); print(head(data$bin_mat.sVar))
      # print("head(data$df.bin.sVar)"); print(head(data$df.bin.sVar))

      super$fit(data) # call the parent class fit method

      # OLD VERSION (when data was still a data.frame)
      # data_t <- data.frame(cbind(self$binirize(data[, self$outvar]), data)) # need a consistent appraoch to creating / handling the dataset
      # super$fit(data_t)
      invisible(self)
    },

    predict = function(newdata) { # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
      if (missing(newdata)) { # ... Do nothing. Predictions for fit data are already saved ...
        return(invisible(self))
      }
      print("predict in continuous summary...")
      # data is of class datnet.sA
      # data$datnetW$dat.sVar -> returns the matrix of sW
      # data$datnetA$dat.sVar -> returns a matrix of sA
      # data$get.df.sW.sA -> returns a data.frame of sWsA

      # NEW VERSION. mat_bin doesn't need to be saved (even though its invisibly returned). 
      # mat_bin is automatically saved in datnet.sW.sA - potentially dangerous!!!
      mat_bin <- newdata$binirize.sVar(self$outvar) # make.bins_mtx_1(make.ordinal(x, intervals), nbins, bin.nms)
      # mat_bin <- newdata$datnetA$binirize.sVar(self$outvar) # make.bins_mtx_1(make.ordinal(x, intervals), nbins, bin.nms)
      print("mat_bin"); print(head(mat_bin))
      super$predict(newdata)

      # OLD VERSION (when data was still a data.frame)
      # newdata_t <- data.frame(cbind(self$binirize(newdata[, self$outvar]), newdata)) # need a consistent appraoch to creating / handling the dataset
      # super$fit(data_t)
      invisible(self)
    },

    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(obsdat.sA) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's)
      assert_that(is.vector(obsdat.sA)) # check A: obsdat.sA is a vector (continuous sA[j])
      print("predictAeqa in continuous summary...")
      
      ordX <- make.ordinal(x = obsdat.sA, intervals = self$bin_intrvls)
      # print("ordX"); print(head(ordX))
      mat_bin <- make.bins_mtx_1(ordX, nbins = self$nbins, bin.nms = self$bin_nms)

      # old way is still possible, since the original datnetA objct was saved as a field:
      # cumprodAeqa <- super$predictAeqa(self$binirize(obsdat.sA))
      cumprodAeqa <- super$predictAeqa(mat_bin)
      invisible(cumprodAeqa)
      # invisible(super$predictAeqa(self$binirize(obsdat.sA)))  # one line alternative to the above
    }
  ),
  active = list(
    cats = function() {seq_len(self$nbins)}
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
