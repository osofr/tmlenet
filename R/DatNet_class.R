#-----------------------------------------------------------------------------
# TO DO: 
#-----------------------------------------------------------------------------
# self$nodes is really unnecessary, all we need to keep are the names of nFnode, Ynode, Anode (for gstar only)
# Defining "equal.mass" bins based on sA under g0 for estimating h_g0 and based on sA under gstar for estimating h_gstar 
  # results in an error, since bins defined under gstar might not include the entire range of sA under g0. 
  # Need a way to fix this so that the intervals are automatically expanded 
  # E.g., add one very very large positive and one very very abs(large) neg value to each interval?)


#-----------------------------------------------------------------------------
# MAJOR CHANGES TO DatNet.sWsA CLASS STRUCTURE (06/18/2015):
#-----------------------------------------------------------------------------
# * gen_sWsA_dat() renamed to make.dat.sWsA() and moved inside DatNet.sWsA class;
# * The observed data is now always stored as $datnetW, $datnetA fields in every DatNet.sWsA object (under g_0 or g_star);
# * DatNet.sWsA$new(datnetW,datnetW) only saves obs data objects datnetW and datnetW as fields;
# * DatNet.sWsA$new() has to be called only twice: once for g_0 (g_N) and once for g_star;
# * $make.dat.sWsA() must be called to create $dat.sWsA ($mat.sVar) - a df / mat of COMBINED sWsA;
# * All binning / interval methods MUST BE called on DatNet.sWsA (NOT $datnetA, $datnetW) (inherited from DatNet);
# * Both copies of DatNet.sWsA are storing datnetA/datnetW by reference - same copy;
# * MOST IMPORTANTLY we can now get rid of the MC sim loop for evaling psi_n. Just use already sampled DatNet.gstar dataset and evaluate psi_n only once.
# * NOTE: Changing datnetA/datnetW in one copy of DatNet.sWsA will result them being changed in the other copy of DatNet.sWsA as well.

#-----------------------------------------------------------------------------
# - DatNet.sWsA$binirize: Current implementation will often create fewer new cats than unique(sVar) for categorical sVar.
  # One way to avoid this is to set bin_by_mass = FALSE for sVar categoricals;
  # Another approach is to collapse the intervals to only unique(intrvls) with a warning;
  # Example: X is cat with 7 unique vals (sum_1mAW2_nets), asked for 7 bins, only 5 non-degenerate bins were created;
  # This is fine for wide format, but will not work when pooling across bins. Pooling requires only non-degenerate bins.
  # Will have to find a way to re-define such bins, so that only 4 final bins are created instead of 7 (with 4th bin being degenerate)
  # [1] "freq count for original variable: "
  #   0   1   2   3   4   5   6 
  # 197 292 234 180  72  18   7 
  # [1] "freq count for transformed ord.sVar: "
  # ord.sVar
  #   2   4   6   7 
  # 197 292 234 277
#-----------------------------------------------------------------------------
# - Is there a way to avoid constant repeating of "self$Var <- Var"? A copy method, s.a., self$copy(arg1, arg2, ...) that will do self$arg1 <- arg1, etc...?
# - See Copy/Close Issue for R6: https://github.com/wch/R6/issues/27 & https://github.com/wch/R6/pull/57 (implemented copy)
  # copytoR6fields <- function(R6obj, ...) {
  #   # put all ... in a named list (unevaluated?)
  #   arglist <- ...
  #   argnames <- names(arglist)
  #   R6obj[[argnames[i]]] <- arglist[[i]]
  # }
#-----------------------------------------------------------------------------
# - Replace all node name references (columns) with indices? See ?subset.data.frame:
  # nl <- as.list(seq_along(data.df))
  # names(nl) <- names(data.df)
  # eval(substitute(node), nl, parent.frame()) -> replace parent.frame() with a reasonable calling envir
## ---------------------------------------------------------------------

is.integerish <- function (x) is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))

#' @importFrom simcausal NetIndClass
#' @importFrom stats as.formula glm na.exclude rbinom 
# @importFrom stats reshape rnorm runif setNames terms.formula
# @importFrom utils head str
# @importFrom graphics legend par plot
NULL


## ---------------------------------------------------------------------
# DETECTING VECTOR TYPES
## ---------------------------------------------------------------------
# sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
detect.col.types <- function(sVar_mat){
  assert_that(is.integerish(getopt("maxncats")) && getopt("maxncats") > 1)
  maxncats <- getopt("maxncats")
  sVartypes <- gvars$sVartypes
  as.list(apply(sVar_mat, 2,
                function(vec) {
                  vec_nomiss <- vec[!gvars$misfun(vec)]
                  nvals <- length(unique(vec_nomiss))
                  if (nvals <= 2L) {
                    sVartypes$bin
                  } else if ((nvals <= maxncats) && (is.integerish(vec_nomiss))) {
                    sVartypes$cat
                  } else {
                    sVartypes$cont
                  }
                }))
}

## ---------------------------------------------------------------------
# Normalizing / Defining bin intervals / Converting contin to ordinal / Converting ordinal to bin indicators
## ---------------------------------------------------------------------
normalize <- function(x) {
  if (abs(max(x) - min(x)) > gvars$tolerr) { # Normalize to 0-1 only when x is not constant
    return((x - min(x)) / (max(x) - min(x)))
  } else {  # What is the thing to do when x constant? Set to abs(x), abs(x)/x or 0???
    return(x)
  }
}
normalize_sVar <- function(sVar_vec) {
  nonmiss_idx <- !gvars$misfun(sVar_vec)
  if (sum(nonmiss_idx) > 0) {
    sVar_vec[nonmiss_idx] <- normalize(sVar_vec[nonmiss_idx])
  }
  sVar_vec
}
normalize_matsVar <- function(sVar_mat) apply(sVar_mat, 2, normalize_sVar)

# Define bin cutt-offs for continuous x:
define.intervals <- function(x, nbins, bin_bymass, bin_bydhist, max_nperbin) {
  x <- x[!gvars$misfun(x)]  # remove missing vals
  nvals <- length(unique(x))

  if (is.na(nbins)) nbins <- as.integer(length(x) / max_nperbin)

  if (nvals < nbins) { # if nbins is too high, for ordinal, set nbins to n unique obs and cancel quantile based interval defns
    nbins <- nvals 
    bin_bymass <- FALSE
  }

  if (abs(max(x) - min(x)) > gvars$tolerr) {  # when x is not constant
    if ((bin_bymass) & !is.null(max_nperbin)) {
      if ((length(x) / max_nperbin) > nbins) nbins <- as.integer(length(x) / max_nperbin)
    }
    intvec <- seq.int(from = min(x), to = max(x) + 1, length.out = (nbins + 1)) # interval type 1: bin x by equal length intervals of 0-1
    # intvec <- c(min(x), sort(runif(n = nbins, min = min(x), max = max(x))), max(x) + 1)
  } else {  # when x is constant, force the smallest possible interval to be at least [0,1]
    intvec <- seq.int(from = min(0L, min(x)), to = max(1L, max(x)), length.out = (nbins + 1))
  }

  if (bin_bymass) {
    intvec <- quantile(x = x, probs = normalize(intvec)) # interval type 2: bin x by mass (quantiles of 0-1 intvec as probs)
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  } else if (bin_bydhist) {
    intvec <- dhist(x, plot = FALSE, nbins=nbins)$xbr
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  }

  # adding -Inf & +Inf as leftmost & rightmost cutoff points to make sure all future data points end up in one of the intervals:
  # intvec <- c(-Inf, min(intvec)-0.01, intvec)
  # intvec <- c(min(intvec) - 0.1, intvec)
  intvec <- c(min(intvec)-1000L, intvec, max(intvec)+1000L)
  # intvec <- c(min(intvec) - 9999999L, min(intvec) - 0.1, intvec, max(intvec) + 0.1, max(intvec) + 9999999L)
  return(intvec) # return(list(intbylen = intvec, intbymass = intvecq))
}

# Turn any x into ordinal (1, 2, 3, ..., nbins) for a given interval cutoffs (length(intervals)=nbins+1)
make.ordinal <- function(x, intervals) findInterval(x = x, vec = intervals, rightmost.closed = TRUE)

# Make dummy indicators for ordinal x (sA[j])
make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms) {
  # Approach: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
  n <- length(x.ordinal)
  cats <- 1 : nbins
  dummies_mat <- matrix(1L, nrow = n, ncol = length(cats))
  for(level in cats[-length(cats)]) {
    subset_Bj0 <- x.ordinal > level
    dummies_mat[subset_Bj0, level] <- 0L
    subset_Bjmiss <- x.ordinal < level
    dummies_mat[subset_Bjmiss, level] <- gvars$misval
  }
  dummies_mat[, cats[length(cats)]] <- gvars$misval
  colnames(dummies_mat) <- bin.nms
  dummies_mat
}

# (4) Approach 2A using model.matrix (a bit slower than manual)
# make.bins_mtx_2A = function(x.ordinal, self) { 
#   # Creates B_j that jumps to 1 and then back to 0 (degenerate) and uses category 1 as ref
#   # First need to make a factor from xcat; adds a column of 1's (ref)
#   n <- length(x.ordinal)
#   cats = self$cats
#   x.ordinalf <- factor(x.ordinal);
#   dummies_mat <- model.matrix(~x.ordinalf);
#   dummies_mat
# }
# (4) Approach 2B (manual vs. of 2A, still slower than 1)
# make.bins_mtx_2B = function(x.ordinal, self) { 
#   # Creates B_j that jumps to 1 and then back to 0 (degenerate), excluding reference cat (last)
#   n <- length(x.ordinal)
#   cats = self$cats
#   sA_nms <- self$sA_nms[-length(self$sA_nms)]
#   dummies2 <- matrix(0L, nrow = n, ncol = length(cats)-1)
#   for(level in cats[-length(cats)]) {
#     dummies2[,level] <- as.integer(x.ordinal == level)
#   }
#   colnames(dummies2) <- sA_nms
#   dummies2
# }

## ---------------------------------------------------------------------
# Class for managing and generating network/summary matrices netVar/sVar where Var is Wnodes or Anode
## ---------------------------------------------------------------------
#' @title Class for storing netW (or netA) (network variables for all W (or A) variables)
#' @docType class
#' @format An R6 class object.
#' @name DatNet
#' @details Following fields are created during initialization
#' \itemize{
#' \item{Kmax} ...
#' \item{netind_cl} ...
#' \item{subset_regs} ...
#' \item{sA_nms} ...
#' \item{sW_nms} ...
#' \item{Kmax} ...
#' }
#' Class for evaluating and storing arbitrary summary measures sVar.
#' The summary measures are evaluated based on the user-specified sVar expressions in sVar.object (sW or sA),
#' in the environment of the input data.frame (Odata). 
#' The evaluated summary measures from sVar.object are stored as a matrix (self$mat.sVar).
#' Contains methods for replacing missing values with default in gvars$misXreplace.
#' Also contains method for detecting /setting sVar variable type (binary, categor, contin).
#' For continous sVar this class provides methods for detecting / setting bin intervals, normalization, disretization and construction of bin indicators.
#' @importFrom assertthat assert_that is.count is.flag
#' @export
DatNet <- R6Class(classname = "DatNet",
  portable = TRUE,
  class = TRUE,
  public = list(
    Kmax = integer(),          # max n of Friends in the network
    nFnode = "nF",
    addnFnode = FALSE,         # Flag to add Fnode to predictors mat / df output
    netind_cl = NULL,          # class NetIndClass object holding $NetInd_k network matrix
    ord.sVar = NULL,           # Ordinal (cat) transform for continous sVar
    active.bin.sVar = NULL,    # name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    Odata = NULL,              # data.frame used for creating the summary measures in mat.sVar, saved each time make.sVar called
    # dat.bin.sVar = NULL,     # (MOVED TO AN ACT BIND) points to self$mat.bin.sVar
    mat.bin.sVar = NULL,       # temp storage mat for bin indicators on currently binarized continous sVar (from self$active.bin.sVar)
    # mat.netVar = NULL,         # NOT DONE mat of network VarNode vals (+ VarNode itself) for each node in VarNodes
    # dat.netVar = NULL,       # (MOVED TO ACT BIND) df of network node values (+ node column itself) for all nodes in VarNodes
    mat.sVar = NULL,           # Matrix storing all evaluated sVars, with named columns
    # dat.sVar = NULL,         # (MOVED TO ACT BIND) Matrix of summary meaure values for each summary measure expression in sVar
    sVar.object = NULL,        # Define_sVar object that contains / evaluates sVar expressions    
    type.sVar = NULL,          # named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    norm.c.sVars = FALSE,      # flag = TRUE if want to normalize continous covariates

    # #todo 37 (DatNet) +0: NEED TO SORT OUT WHEN nOdata is needed....
    # nobs = NA_integer_,      # n of samples in the OBSERVED (original) data
    nOdata = NA_integer_,      # n of samples in the OBSERVED (original) data

    initialize = function(netind_cl, nodes, nFnode, addnFnode = FALSE, ...) {
      assert_that(is.flag(addnFnode))
      self$netind_cl <- netind_cl
      self$Kmax <- netind_cl$Kmax
      if (!missing(nFnode)) self$nFnode <- nFnode
      self$addnFnode <- addnFnode
      if (!missing(nodes)) self$nodes <- nodes
      invisible(self)
    },

    # **********************
    # Define summary measures sVar
    # **********************
    # type.sVar acts as a flag: only detect types when !is.null, addnFnode = TRUE
    make.sVar = function(Odata, sVar.object = NULL, type.sVar = NULL, norm.c.sVars = FALSE) {
      assert_that(is.data.frame(Odata))
      self$nOdata <- nrow(Odata)
      self$Odata <- Odata

      if (is.null(sVar.object)) {
        stop("Not Implemented. To Be replaced with netVar construction when sVar.object is null...")
      }

      self$sVar.object <- sVar.object
      if (self$addnFnode) { nFnode <- self$nFnode } else { nFnode <- NULL }
      self$mat.sVar <- sVar.object$get.mat.sVar(data.df = Odata, netind_cl = self$netind_cl, addnFnode = nFnode)

      # MAKE def_types_sVar an active binding? calling self$def_types_sVar <- type.sVar assigns, calling self$def_types_sVar defines.
      self$def_types_sVar(type.sVar) # Define the type of each sVar[i]: bin, cat or cont

      # normalize continuous and non-missing sVars, overwrite their columns in mat.sVar with normalized [0,1] vals
      if (norm.c.sVars) {
        self$norm.c.sVars <- norm.c.sVars
        self$norm_c_sVars()
      }
      invisible(self)
    },

    # *** MAKE A PRIVATE METHOD ***
    # Define the type (class) of each summary measure: bin, cat or cont
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar)
    # otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    def_types_sVar = function(type.sVar = NULL) {
      # Detect the type of each sVar[i]: gvars$sVartypes$bin,  gvars$sVartypes$cat, gvars$sVartypes$cont
      if (is.null(type.sVar)) {
        self$type.sVar <- detect.col.types(self$dat.sVar)
      } else {
        n.sVar <- length(self$names.sVar)
        len <- length(type.sVar)
        assert_that((len == n.sVar) || (len == 1L))
        if (len == n.sVar) {
          assert_that(is.list(type.sVar))
          assert_that(all(names(type.sVar) %in% self$names.sVar))
        } else {
          assert_that(is.string(type.sVar))
          type.sVar <- as.list(rep(type.sVar, n.sVar))
          names(type.sVar) <- self$names.sVar
        }
        self$type.sVar <- type.sVar
      }
      invisible(self)
    },

    # *** MAKE A PRIVATE METHOD ***
    # Normalize continuous sVars # This could be memory-costly
    norm_c_sVars = function() {
      names.c.sVar <- self$names.c.sVar
      if (length(names.c.sVar) == 0L) return(invisible(self))

      if (self$norm.c.sVars && (length(names.c.sVar) > 0)) {
        for (name.c.sVar in names.c.sVar) {
          self$mat.sVar[, name.c.sVar] <- normalize_sVar(self$mat.sVar[, name.c.sVar])
        }
      }
      invisible(self)
    },

    # #todo 18 (DatNet, DatNet.sWsA) +0: (OPTIONAL) ENABLE ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO DatNet
    # #todo 25 (DatNet, add_det) +0: Need to save det node flags as a separate mat, can't add them to sVar since all sVars 
    # will be automatically added to A ~ predictors
    # add_deterministic = function(Odata, userDETcol) {
    #   determ.g_user <- as.vector(Odata[,userDETcol]) # get deterministic As for the entire network of each unit (set by user)
    #   # determ.gvals_user <- Odata[,AnodeDET] # add values to be assigned to deterministic nodes (already have: netA)
    #   determ_cols_user <- .f.allCovars(k = self$Kmax, NetInd_k = self$netind_cl$NetInd_k, Var = determ.g_user,
    #                                     VarNm = "determ.g_true", misval = gvars$misval)
    #   # determ_cols <- (determ_cols_user | determ_cols_Friend)
    #   determ_cols <- determ_cols_user
    #   # THIS IS A VERY BAD WAY TO DO THAT, REMOVING
    #   # self$mat.sVar <- cbind(determ_cols, self$mat.sVar)
    #   determ_cols_type <- as.list(rep.int(gvars$sVartypes$bin, ncol(determ_cols)))
    #   names(determ_cols_type) <- colnames(determ_cols)
    #   self$type.sVar <- c(determ_cols_type, self$type.sVar)
    #   invisible(self)
    # },

    # # Replaces all misval values in mat.netVar / mat.sVar with gvars$misXreplace
    # fixmiss_netVar = function() {
    #   self$mat.netVar[gvars$misfun(self$mat.netVar)] <- gvars$misXreplace
    #   invisible(self)
    # },
    fixmiss_sVar = function() {
      self$mat.sVar[gvars$misfun(self$mat.sVar)] <- gvars$misXreplace
      invisible(self)
    },

    # --------------------------------------------------
    # Methods for directly handling one continous/categorical sVar in self$mat.sVar;
    # No checking of incorrect input is performed, use at your own risk!
    # --------------------------------------------------
    norm.sVar = function(name.sVar) { normalize_sVar(self$dat.sVar[, name.sVar]) },  # return normalized 0-1 sVar
    set.sVar = function(name.sVar, new.sVar) { self$mat.sVar[, name.sVar] <- new.sVar },
    get.sVar = function(name.sVar) { self$dat.sVar[, name.sVar] },
    set.sVar.type = function(name.sVar, new.type) { self$type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { self$type.sVar } else { self$type.sVar[[name.sVar]] } }

  ),

  active = list(
    # names.netVar = function() { colnames(self$dat.netVar) },
    names.sVar = function() { colnames(self$dat.sVar) },
    names.c.sVar = function() { names(self$type.sVar[self$type.sVar %in% gvars$sVartypes$cont]) }, # names of cont sVars
    ncols.netVar = function() { length(self$names.netVar) },
    ncols.sVar = function() { length(self$names.sVar) },

    # dat.netVar = function() { self$mat.netVar },
    dat.sVar = function() { self$mat.sVar },
    dat.bin.sVar = function() { self$mat.bin.sVar },

    # emptydat.netVar = function() {self$mat.netVar <- NULL },    # wipe out mat.netVar
    emptydat.sVar = function() { self$mat.sVar <- NULL },         # wipe out mat.sVar
    emptydat.bin.sVar = function() { self$mat.bin.sVar <- NULL }, # wipe out binirized mat.sVar
    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        private$.nodes <- nodes
      }
    }
  ),
  private = list(
    .nodes = list()           # names of the nodes in the data (Anode, Ynode, nFnode, etc..)
  )
)

## ---------------------------------------------------------------------
# Class for managing and generating all the summary datasets sWsA used in fitting \bar{h}^* or \bar{h}_0
# DatNet.sWsA is the only way to access data in the entire package 
# This class gets passed on to SummariesModel functions: $fit(), $predict() and $predictAeqa()
# #todo 33 (DatNet.sWsA) +0: rename datnetW, datnetA to O.datnetW, O.datnetA for clarity
## ---------------------------------------------------------------------
#' @title Class for managing combined summary measures sW & sA from DatNet class.
#' @docType class
#' @format An R6 class object.
#' @name DatNet.sWsA
#' @details Following fields are created during initialization
#' \itemize{
#' \item{nodes} ...
#' \item{datnetW} ...
#' \item{datnetA} ...
#' \item{YnodeVals} ...
#' \item{det.Y} ...
#' \item{Kmax} ...
#' }
#' Methods for combining, subsetting, discretizing & binirizing of summary measures in sW & sA. Inherits from DatNet class. 
#' The combined dataset of all (sW, sA) summary measures returned by DatNet is stored as a matrix under self$mat.sVar
#' @importFrom assertthat assert_that is.count is.flag
#' @export
DatNet.sWsA <- R6Class(classname = "DatNet.sWsA",
  inherit = DatNet,
  portable = TRUE,
  class = TRUE,
  public = list(
    datnetW = NULL,   # *** RENAME TO O.datnetW for clarity ***
    datnetA = NULL,   # *** RENAME TO O.datnetA for clarity ***
    YnodeVals = NULL, # Values of the binary outcome (Ynode) in observed data where det.Y = TRUE obs are set to NA
    det.Y = NULL,     # Logical vector, where YnodeVals[det.Y==TRUE] are deterministic (0 or 1)
    p = 1,
    # **********
    # dat.sVar - (inherited act bind): now points to combine mat.sVar of above cbind(dat.sW, dat.sA)
    # dat.bin.sVar - (inherited act bind): points to mat of binned sVar (currently selected) which is mat.bin.sVar
    # active.bin.sVar - name of the currently selected binirized sVar
    # this keeps ALL methods and active bindings of DatNet valid in DatNet.sWsA for this combined data mat
    # **********

    initialize = function(datnetW, datnetA, YnodeVals, det.Y, ...) {
      assert_that("DatNet" %in% class(datnetW))
      assert_that("DatNet" %in% class(datnetA))
      self$datnetW <- datnetW
      self$datnetA <- datnetA
      self$netind_cl <- datnetW$netind_cl
      self$Kmax <- self$netind_cl$Kmax
      # re-assign nodes object if it already exists in datnetW
      if (length(datnetW$nodes) > 0) self$nodes <- datnetW$nodes
      # self$nodes <- datnetW$nodes
      if (!missing(YnodeVals)) self$addYnode(YnodeVals = YnodeVals, det.Y = det.Y)
      invisible(self)
    },

    addYnode = function(YnodeVals, det.Y) {
        if (missing(det.Y)) det.Y <- rep.int(FALSE, length(YnodeVals))
        self$noNA.Ynodevals <- YnodeVals  # Adding actual observed Y as protected (without NAs)
        self$YnodeVals <- YnodeVals
        self$YnodeVals[det.Y] <- NA       # Adding public YnodeVals & setting det.Y values to NA
        self$det.Y <- det.Y
    },
 
    evalsubst = function(subsetexpr, subsetvars) { # Eval the expression (in the environment of the data.frame "data" + global constants "gvars"):      
      # Could also do evaluation in a special env with a custom subsetting fun '[' that will dynamically find the correct dataset that contains 
      # sVar.name (dat.sVar or dat.bin.sVar) and will return sVar vector
      if (missing(subsetexpr)) {
        assert_that(!missing(subsetvars))
        res <- rep.int(TRUE, self$nobs)

        for (subsetvar in subsetvars) {
          # print("evaluating subset vector for var: " %+% subsetvar)
          # *) find the var of interest (in self$dat.sWsA or self$dat.bin.sVar), give error if not found
          sVar.vec <- self$get.outvar(var = subsetvar)
          assert_that(!is.null(sVar.vec))
          # *) reconstruct correct expression that tests for missing values
          res <- res & (!gvars$misfun(sVar.vec))
        }
        return(res)

      } else {
        if (is.logical(subsetexpr)) {
          return(subsetexpr)
        } else {
          # ******************************************************
          # THIS WAS A BOTTLENECK
          # for 500K w/ 1000 bins: 4-5sec
          # ******************************************************
          # REPLACING WITH env that is made of data.frames instead of matrices
            eval.env <- c(data.frame(self$dat.sWsA), data.frame(self$dat.bin.sVar), as.list(gvars))
            res <- try(eval(subsetexpr, envir = eval.env, enclos = baseenv())) # to evaluate vars not found in data in baseenv()
          return(res)
        }
      }
    },

    get.dat.sWsA = function(rowsubset = TRUE, covars) { # return a data.frame with covars (predictors)
      dat.bin.sVar <- self$dat.bin.sVar
      sel.sWsA <- TRUE
      sel.binsA = NULL # columns to select from binned continuos var matrix (if it was previously constructed)
      if (!missing(covars)) {
        sel.sWsA <- colnames(self$dat.sWsA)[(colnames(self$dat.sWsA) %in% covars)]
        if (!is.null(dat.bin.sVar)) {
          sel.binsA <- colnames(dat.bin.sVar)[(colnames(dat.bin.sVar) %in% covars)]
        }
      }
      dfsel <- self$dat.sWsA[rowsubset, sel.sWsA, drop = FALSE]
      if (length(sel.binsA)>0) {
        dfsel <- cbind(dfsel, dat.bin.sVar[rowsubset, sel.binsA, drop = FALSE])
      }
      return(dfsel)
    },

    get.outvar = function(rowsubset = TRUE, var) {
      if (length(self$nodes) < 1) stop("DatNet.sWsA$nodes list is empty!")
      if (var %in% self$names.sWsA) {
        self$dat.sWsA[rowsubset, var]
      } else if (var %in% colnames(self$dat.bin.sVar)) {
        self$dat.bin.sVar[rowsubset, var]
      } else if ((var %in% self$nodes$Ynode) && !is.null(self$YnodeVals)) {
        self$YnodeVals[rowsubset]
      } else {
        # print("var: "); print(var)
        # print("self$YnodeVals: "); print(head(self$YnodeVals))
        stop("requested variable " %+% var %+% " does not exist in DatNet.sWsA!")
        # NULL
      }
    },

    copy.sVar.types = function() {
      self$type.sVar <- c(self$datnetW$type.sVar, self$datnetA$type.sVar)
    },

    # ------------------------------------------------------------------------------------------------------------
    # MOVED THESE METHODS TO DatNet.sWsA class
    # ------------------------------------------------------------------------------------------------------------
    # Need to find a way to over-ride nbins for categorical vars (allowing it to be set to more than gvars$maxncats)!
    # Return names of bin indicators for sVar:
    bin.nms.sVar = function(name.sVar, nbins) { name.sVar%+%"_"%+%"B."%+%(1L:nbins) }, 
    pooled.bin.nm.sVar = function(name.sVar) { name.sVar %+% "_allB.j" },
    # #todo 71 (DatNet.sWsA) +0: *** NOTE *** When sVar is cat might be better to set bin_bymass = FALSE to avoid collapsing of categories for sVar
    detect.sVar.intrvls = function(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin) {
      int <- define.intervals(x = self$dat.sVar[, name.sVar], nbins = nbins, bin_bymass = bin_bymass, bin_bydhist = bin_bydhist, max_nperbin = max_nperbin)
      if (length(unique(int)) < length(int)) {
        message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+% 
                (length(int)-1) %+% " to " %+% (length(unique(int))-1) %+% " due to too few obs.")
        print("old intervals: "); print(int)
        int <- unique(int)
        print("new intervals: "); print(int)
      }
      return(int)
    },
    # create a vector of ordinal (categorical) vars out of cont. sVar vector:
    discretize.sVar = function(name.sVar, intervals) {
      self$ord.sVar <- make.ordinal(x = self$dat.sVar[, name.sVar], intervals = intervals)
      invisible(self$ord.sVar)
    },
    # return matrix of bin indicators for ordinal sVar:
    binirize.sVar = function(name.sVar, intervals, nbins, bin.nms) {
      self$active.bin.sVar <- name.sVar
      self$mat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$discretize.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms)
      invisible(self$mat.bin.sVar)
    },
    # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bw = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      intrvls.width <- diff(intervals)
      intrvls.width[intrvls.width <= gvars$tolerr] <- 1
      ord.sVar_bw <- intrvls.width[self$ord.sVar]
      return(ord.sVar_bw)
    },
   # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bwdiff = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      # intrvls.width <- diff(intervals)
      # intrvls.width[intrvls.width <= gvars$tolerr] <- 1
      ord.sVar_leftint <- intervals[self$ord.sVar]
      diff_bw <- self$dat.sVar[, name.sVar] - ord.sVar_leftint
      return(diff_bw)
    },

    # This function returns mat.sVar, which is a matrix that combines all sW and sA summary measures;
    # Odata is only needed for evaluating new sA (!is.null(f.g_name));
    # When !is.null(f.g_name) create p new datnetA.gstar's (n obs at a time), which are not saved separately (only combined);
    # When is.null(f.g_name), returns combined cbind(sW, sA) for observed O.datnetW, O.datnetA;
    # TO ADD: Consider passing ahead a total number of sA that will be created by DatNet class (need it to pre-allocate self$dat.sWsA);
    make.dat.sWsA = function(p = 1, f.g_name = NULL, f.g_args = NULL, sA.object = NULL)  {
      datnetW <- self$datnetW
      datnetA <- self$datnetA
      assert_that(is.count(p))
      self$p <- p
      nobs <- datnetW$nOdata
      # Copy variable detected types (bin, cat or contin) from the observed data classes (datnetW, datnetA) to self:
      self$copy.sVar.types()
      if (is.null(f.g_name)) {  # set df.sWsA to observed data (sW,sA) if g.fun is.null
        df.sWsA <- cbind(datnetW$dat.sVar, datnetA$dat.sVar) # assigning summary measures as data.frames:
      } else {  # need to sample A under f.g_name (gstar or known g0), possibly re-evaluate sW from O.datnetW
        if (is.null(self$nodes$Anode)) stop("Anode was not appropriately specified and is null; can't replace observed Anode with that sampled under g_star")
        Odata <- datnetW$Odata
        # Will not be saving this object datnetA.gstar as self$datnetA (i.e., keeping a old pointer to O.datnetA)
        datnetA.gstar <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes)
        # datnetA.gstar <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes, VarNodes = self$nodes$Anode)
        df.sWsA <- matrix(nrow = (nobs * p), ncol = (datnetW$ncols.sVar + datnetA$ncols.sVar))  # pre-allocate result matx sWsA
        colnames(df.sWsA) <- self$names.sWsA
        for (i in seq_len(p)) {
          # *** f.g_name can only depend on covariates in datnetW$dat.sVar ***
          # if Anode is continuous, just call f.gen.probA.star:
          A.gstar <- f.gen.A.star.cont(k = self$Kmax, df_AllW = cbind(datnetW$dat.sVar,datnetA$dat.sVar), fcn_name = f.g_name, f_args = f.g_args)
          Odata[, self$nodes$Anode] <- A.gstar # replace A under g0 in Odata with A^* under g.star:
          datnetA.gstar$make.sVar(Odata = Odata, sVar.object = sA.object) # create new summary measures sA (under g.star)
          # Assigning the summary measures to one output data matrix:
          df.sWsA[((i - 1) * nobs + 1):(nobs * i), ] <- cbind(datnetW$dat.sVar, datnetA.gstar$dat.sVar)[, ]
        }
      }
      self$mat.sVar <- df.sWsA
      invisible(self)
    }
  ),
  active = list(
    dat.sWsA = function() { self$mat.sVar }, # NO LONGER NEEDED, KEPT FOR COMPATIBILITY
    names.sWsA = function() { c(self$datnetW$names.sVar, self$datnetA$names.sVar) },
    nobs = function() { nrow(self$dat.sWsA) },
    noNA.Ynodevals = function(noNA.Yvals) {
      if (missing(noNA.Yvals)) return(private$protected.YnodeVals)
      else private$protected.YnodeVals <- noNA.Yvals
    },
    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        if (length(self$nodes)>0) message("warning: overwriting non-empty self$nodes in DatNet.sWsA")
        private$.nodes <- nodes
        # propagate new nodes to parent objects:
        if (length(self$datnetW$nodes) < 1) self$datnetW$nodes <- nodes
        if (length(self$datnetA$nodes) < 1) self$datnetA$nodes <- nodes
      }
    }    
  ),
  private = list(
    protected.YnodeVals = NULL  # Actual observed values of the binary outcome (Ynode), along with deterministic vals
  )
)