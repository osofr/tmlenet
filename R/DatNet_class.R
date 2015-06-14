library(assertthat)
library(speedglm)
library(R6)
`%+%` <- function(a, b) paste0(a, b)  # cat function
#-----------------------------------------------------------------------------
# TO DO: 
#-----------------------------------------------------------------------------
# - Is there a way to avoid constant repeating of "self$Var <- Var"? A copy method, s.a., self$copy(arg1, arg2, ...) that will do self$arg1 <- arg1, etc...?
# - See Copy/Close Issue for R6: https://github.com/wch/R6/issues/27 & https://github.com/wch/R6/pull/57 (implemented copy)
copytoR6fields <- function(R6obj, ...) {
  # put all ... in a named list (unevaluated?)
  arglist <- ...
  argnames <- names(arglist)
  R6obj[[argnames[i]]] <- arglist[[i]]
}
# - Allow different summary measures sA and sW for g_0 and g_star (h_0 and h_star) 
  # x) Create a class out of fit.hbars()?
# - Automatically determine when values sA[j], sW[j] are degenerate/deterministic:  
  # x) excluded from h regs when sA[j] is an outcome
  # x) replaced with 0 vals when sA[j]/sW[j] are covariates
  # x) enable assignment of probA1 & probA.eq.a for deterministic nodes (including user def. deterministic Anode)
# - Sort out passing/storing (1) data_mtx, (2) new_mtx, (3) X_mat without having to make extra copies
# - Continuous summary measures?

## ---------------------------------------------------------------------
# Very simple class to hold and create NetInd_k, the matrix of network connection indices in observed data, dim = (nobs x Kmax)
## ---------------------------------------------------------------------
NetIndClass <- R6Class("NetIndClass",
  class = FALSE,
  portable = FALSE,
  public = list(
    NetInd_k = NULL,
    nobs = NA_integer_,
    Kmax = NA_integer_,
    initialize = function(Odata, Kmax) {
      nobs <<- nrow(Odata)
      Kmax <<- Kmax
      makeNetInd(Odata = Odata)
    },
    makeNetInd = function(Odata) {
      NetInd_k <<- matrix(0L, nrow = nobs, ncol = Kmax)
    }
  ),
  active = list(
    getNetInd = function() NetInd_k
  )
)
# testdf <- data.frame(a = rnorm(100), b = rnorm(100))
# nettest <- NetIndClass$new(Odata = testdf, Kmax = 5)
# nettest$getNetInd
# nettest$getNetInd <- NULL

## ---------------------------------------------------------------------
# Class for managing and generating network/summary matrices netVar/sVar where Var is Wnodes or Anode
## ---------------------------------------------------------------------
# TO DO: 
# x) See if new() can be made to sample Odata[,Varnode] if f.g_name is supplied... 
  # This would allow making DatNet fully closed in, that is, its the only class where netVar and sVar data are contructed, given Odata and f.g_name...
# x) See if some of new() args can be dropped: (Odata, NetInd_k, Kmax, nodes, VarNodes, AddnFnode = FALSE, misValRepl = FALSE, ...) {

#' @title Class for storing netW (or netA) (network variables for all W (or A) variables)
#' @docType class
#' @format An R6 class object.
#' @name DatNetW
#' @details Following fields are created during initialization
#' \itemize{
#' \item{nodes} ...
#' \item{subset_regs} ...
#' \item{sA_nms} ...
#' \item{sW_nms} ...
#' \item{Kmax} ...
#' }
#' Class for constructing and storing one network dataset netW for W or netA for A, as well as the summary data set sW or sA.
##' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
##' @export
DatNet <- R6Class(classname = "DatNet",
  portable = TRUE,
  class = TRUE,
  public = list(
    Kmax = integer(),          # max n of Friends in the network
    nodes = list(),            # names of the nodes in the data (Anode, Wnodes, nFnode, etc..)
    NetInd_k = NULL,           # class NetIndClass object holding NetInd_k network matrix
    VarNodes = NULL,
    AddnFnode = FALSE,
    misValRepl = FALSE,

    dat.netVar = NULL,         # matrix of network node values (+ node column itself) for all nodes in VarNodes
    names.netVar = character(),# column names of dat.netVar. Consider turning it into an active field (with arg)

    dat.sVar = NULL,           # matrix of summary meaure values for each summary measure expression in sVar
    names.sVar = character(),  # summary meaure names. Consider turning it into an active field (with arg)
    nobs = NA_integer_,        # n of samples in the OBSERVED (original) data

    initialize = function(Odata, NetInd_k, Kmax, nodes, VarNodes, AddnFnode = FALSE, misValRepl = FALSE, ...) {
      # Ovector,
      assert_that(!is.null(VarNodes)); assert_that(is.list(VarNodes)||is.character(VarNodes))
      self$VarNodes <- VarNodes # self$VarNodes <- nodes$Wnodes # or: # self$VarNodes <- nodes$Anode
      # a better alternative is: # self$VarNodes <- nodes[[VarNodeType]] # VarNodeType <- "Wnodes"; VarNodeType <- "Anode"

      assert_that(is.count(Kmax))
      assert_that(is.flag(AddnFnode))
      assert_that(is.flag(misValRepl))
      
      self$NetInd_k <- NetInd_k   # include or not?
      self$Kmax <- Kmax
      self$nodes <- nodes
      # if (!missing(Odata)) {
        assert_that(is.data.frame(Odata)) 
        self$nobs <- nrow(Odata)
      # } else {
      #   assert_that(!missing(Ovector))
      #   assert_that(is.vector(Ovector))
      #   assert_that(!AddnFnode)
      #   assert_that(length(VarNodes)==1L)
      #   self$nobs <- length(Ovector)
      # }

      self$AddnFnode <- AddnFnode
      self$misValRepl <- misValRepl

      self$names.netVar <- netvar2(unlist(self$VarNodes), (0L:self$Kmax))
      if (self$AddnFnode) self$names.netVar <- c(self$names.netVar, self$nodes$nFnode)

      # Build network vectors: (W, W_netF_1, ..., W_netF_k) for each W in Wnodes by PRE-ALLOCATING netW_full:      
      self$dat.netVar <- matrix(0L, nrow = self$nobs, ncol = self$ncols.netVar)
      colsperVar <- (self$Kmax + 1)
      for (idxVar in seq(self$VarNodes)) {
        Varnode <- self$VarNodes[[idxVar]]
        self$dat.netVar[, ((idxVar - 1) * colsperVar + 1) : (idxVar * colsperVar)] <- 
              .f.allCovars(k = self$Kmax, NetInd_k = self$NetInd_k, Var = Odata[,Varnode], VarNm = Varnode, misval = gvars$misval)
      }
      if (self$AddnFnode) self$dat.netVar[, self$ncols.netVar] <- as.matrix(Odata[, self$nodes$nFnode, drop=FALSE])
      colnames(self$dat.netVar) <- self$names.netVar

      if (self$misValRepl) self$fixmiss_netVar() # replace missing vals with gvars$misXreplace

      invisible(self)
    },

    # II) NOT FULLY IMPLEMENTED: For now can only specify summary measure as subsets of netVar columns
    # APPLY THE SUMMARY MEASURE FUNCTIONS / EXPRESSION TO netW_full to OBTAIN sW columns
    # SELECT ONLY sW columns in hform_g0 and hfrom_gstar or use all?
    make_sVar = function(names.sVar = NULL) {
      if (is.null(names.sVar)) {
        names.sVar <- self$names.netVar
      }
      assert_that(is.character(names.sVar))
      assert_that(all(names.sVar %in% self$names.netVar))
      self$names.sVar <- names.sVar
      self$dat.sVar <- self$dat.netVar[, self$names.sVar]
      invisible(self)
    },

    # (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sW:
    add_deterministic = function(Odata, userDETcol) {
      determ.g_user <- as.vector(Odata[,userDETcol]) # get deterministic As for the entire network of each unit (set by user)
      # determ.gvals_user <- Odata[,AnodeDET] # add values to be assigned to deterministic nodes (already have: netA)
      determ_cols_user <- .f.allCovars(k = self$Kmax, NetInd_k = self$NetInd_k, Var = determ.g_user, 
                                        VarNm = "determ.g_true", misval = gvars$misval)
      # determ_cols <- (determ_cols_user | determ_cols_Friend)
      determ_cols <- determ_cols_user
      self$dat.sVar <- cbind(determ_cols, self$dat.sVar)
      self$names.sVar <- c(colnames(determ_cols), self$names.sVar)
      invisible(self)
    },
    # REPLACE ALL misval values in sW with gvars$misXreplace (OR DO IT ON THE ACTUAL DATASET WHEN SUBSETTING IS PERFORMED)
    fixmiss_netVar = function() {
      self$dat.netVar[gvars$misfun(self$dat.netVar)] <- gvars$misXreplace
      invisible(self)
    },
    fixmiss_sVar = function() {
      self$dat.sVar[gvars$misfun(self$dat.sVar)] <- gvars$misXreplace
      invisible(self)
    }
  ),
  active = list(
    ncols.netVar = function() { length(self$names.netVar) },
    ncols.sVar = function() { length(self$names.sVar) },
    emptydat.netVar = function() {self$dat.netVar <- NULL},
    emptydat.sVar = function() { self$dat.sVar <- NULL }
  ),
  private = list(
    # dat.netVar = NULL,
    # dat.sVar = NULL,
    placeholder = list()
  )
)


## ---------------------------------------------------------------------
# Class for managing and generating all the summary datasets sWsA for fitting \bar{h}^* or \bar{h}_0
## ---------------------------------------------------------------------
# TO DO: 
# x) MAKE DatSummaries inherit from DatNet and DO NOT REPLICATE ANY FIELDS THAT ALREADY EXIST IN DatNet!

# #' @title Class for storing summary measure datasets.
# #' @docType class
# #' @format An R6 class object.
# #' @name DatSummaries
# #' @details Following fields are created during initialization
# #' \itemize{
# #' \item{nodes} ...
# #' \item{subset_regs} ...
# #' \item{sA_nms} ...
# #' \item{sW_nms} ...
# #' \item{Kmax} ...
# #' }
# #' Class for constructing and storing the summary measure dataset (netVar, sVar), where Var is either A or W and sVar is either sW or sA.
# ##' @importFrom R6 R6Class
# #' @importFrom assertthat assert_that
# ##' @export
# DatSummaries <- R6Class(classname = "DatSummaries",
#   # inherit = DatNet,
#   portable = TRUE,
#   class = TRUE,
#   public = list(
#     Kmax = integer(),       # max n of Friends in the network
#     nodes = list(),    # names of the nodes in the data (Anode, Wnodes, nFnode, etc..)
#     f.g_name = character(),
#     f.g_args = list(),
#     Aobs.vec = NA_integer_, # vector of observed A values
#     NetInd_k = NULL,
#     dat.W = NULL,           # an object of class DatNet that stores baseline netW and sW matrices
#     # dat.netA = NULL, # dat.sW = NULL, # dat.sA = NULL,
#     dat.sWsA = NULL,
#     len.sW = NA_integer_,    # number of summary measures used for sW
#     len.sA = NA_integer_,    # number of summary measures used for sW
#     names.sW = character(),  # all sW names
#     names.sA = character(),  # all sA names
#     p = NA_integer_,        # number of repeat samples from g, for each sample of size nrow(data)
#     nobs = NA_integer_,        # n of samples in the OBSERVED (original) data
    
#     initialize = function(Odata, Kmax, nodes, NetInd_k, f.g_name, f.g_args, dat.W, ...) {
#       self$f.g_name <- f.g_name
#       self$f.g_args <- f.g_args
#       assert_that(is.count(Kmax))
#       self$Kmax <- Kmax
#       # pass this directly to gen_sWsA_dat fun:
#       self$nodes <- nodes
#       self$Aobs.vec <- Odata[, nodes$Anode]
#       self$NetInd_k <- NetInd_k # include or not?
#       assert_that("dat.W" %in% class(dat.W))
#       self$dat.W <- dat.W
#       # assert_that(is.character(sA_nms)); assert_that(is.character(sW_nms))
#       # self$sA_nms <- sA_nms; self$sW_nms <- sW_nms;
#       self$nobs <- nrow(Odata)
#       invisible(self)
#     },

#     #---------------------------------------------------------------------------------
#     # SAMPLE A LARGE DATASET of cY's for given functions g* or g0, of size p*nobs for some p
#     # Returns observed network data if is.null(f.g_name)
#     # NEED TO ADD: pass ahead a total number of sA that will be created by DatNet class (need it to pre-allocate self$dat.sWsA)
#     # Current structure requires building sA twice, once for observed data and once for g_0 when g_0 unknown. This can be expensive. 
#     #---------------------------------------------------------------------------------
#     gen_sWsA_dat = function(p = 1)  {
#      # ... see fit.h()
#       invisible(self)
#     }
#   ),
#   active = list(
#     getAvec = function() {  # either sample A vector from f.g_name or return observed A vector
#       if (is.null(self$f.g_name)) {
#         self$Aobs.vec
#       } else {
#         f.gen.A.star(self$Kmax, self$dat.W$dat.netVar, self$f.g_name, self$f.g_args)
#       }
#       # dat.sA = function(){
#       #   self$dat.sWsA[, self$names.sA]
#       # }
#     },
#     emptynetdata = function() {
#       # private$Y_vals <- NULL
#       # private$X_mat <- NULL
#       # self$subset <- NULL # wipe out subset as well or leave it alone?
#     }
#   ),

#   private = list(
#     placeholder = list()
#   )
# )
