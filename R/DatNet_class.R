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
# Changes to Subsetting / Delayed eval (DatNet.sWsA$evalsubst)
#-----------------------------------------------------------------------------
# x) Delayed eval. is inside DatNet.sWsA
# x) Delayed evaluation of subsets implies that data / newdata have to be data.frames at the time of subseting.
# x) **** Since data is now stored as matrix that means 3 conversions matrix -> df for each call to $evalsubst **** 
# x) See implementation of subset.data.frame and subset.matrix

#-----------------------------------------------------------------------------
# TO DO: 
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

# - DatNet.sWsA$evalsubst: Instead of expressions need to call pre-defined subsetting functions?
# - DatNet.sWsA$evalsubst: Could try the same approach as in Define_sVar, using custom '[' to select rows and cols based on subset_expr

# - Is there a way to avoid constant repeating of "self$Var <- Var"? A copy method, s.a., self$copy(arg1, arg2, ...) that will do self$arg1 <- arg1, etc...?
# - See Copy/Close Issue for R6: https://github.com/wch/R6/issues/27 & https://github.com/wch/R6/pull/57 (implemented copy)
  # copytoR6fields <- function(R6obj, ...) {
  #   # put all ... in a named list (unevaluated?)
  #   arglist <- ...
  #   argnames <- names(arglist)
  #   R6obj[[argnames[i]]] <- arglist[[i]]
  # }
# - Replace all node name references (columns) with indices? See ?subset.data.frame:
  # nl <- as.list(seq_along(data.df))
  # names(nl) <- names(data.df)
  # eval(substitute(node), nl, parent.frame()) -> replace parent.frame() with a reasonable calling envir
## ---------------------------------------------------------------------


## ---------------------------------------------------------------------
# Class holds and creates NetInd_k, the matrix of network connection indices in Odata of dim = (nobs x Kmax)
# Also calculates a vector nF - number of friends for each unit
# When class = FALSE, a pointer "self" is still created, but parent.env(self) is the enclosing environment
## ---------------------------------------------------------------------
NetIndClass <- R6Class("NetIndClass",
  class = FALSE,
  portable = FALSE,
  public = list(
    NetInd_k = matrix(),       # matrix (n x Kmax) of network (friend) indices (rows) in observed data
    # NetVec_l = list(),       # (NOT USED) list (n) of network (friend) indices (rows) in observed data
    nF = integer(),            # number of friends, integer vector of length n
    nobs = NA_integer_,        # n observations
    Kmax = NA_integer_,        # max number of friends
    IDnode = NULL,             # name of the column in Odata with unit ids
    NETIDnode = character(),   # name of the column in Odata with a list of friend ids (as a string, with ids separated by sep)
    sep = ' ',                 # character separating two friend ids from column NETIDnode

    initialize = function(Odata, Kmax, IDnode = NULL, NETIDnode, sep = ' ') {
      assert_that(is.data.frame(Odata))
      nobs <<- nrow(Odata)

      assert_that(is.count(Kmax))
      Kmax <<- Kmax

      if (!is.null(IDnode)) {
        assert_that(is.character(IDnode))
        IDnode <<- IDnode
      }

      assert_that(is.character(NETIDnode))
      NETIDnode <<- NETIDnode

      assert_that(is.string(sep))
      sep <<- sep

      # could pre-allocate NetInd_k mat: # NetInd_k <<- matrix(0L, nrow = nobs, ncol = Kmax)
      resNets <- getNetInd(data = Odata, Kmax = Kmax, IDnode = IDnode, NETIDnode = NETIDnode, sep = sep)
      
      nF <<- resNets$nF
      NetInd_k <<- resNets$NetInd_k
      invisible(self)
    },

    assignNetInd = function(NetInd_k) {
      NetInd_k[,] <<- NetInd_k
    },

    #------------------------------------------------------------------------------
    # Netwk ids strings to list of friend indices (NetVec) and matrix of columns of friend indices (NetInd_k)
    #------------------------------------------------------------------------------
    # data - input data.frame with a column for friend IDs
    # Kmax - max number of friends
    # IDnode - column name in data for subject IDs
    # NETIDnode - name of the node in data that contains friend IDs as a character string
    # sep - character symbol separating two friend IDs in data[i, NETIDnode] for observation i
    getNetInd = function(data, Kmax, IDnode = NULL, NETIDnode, sep = ' ') {

      # Turn string of IDs into a vector, trim extra spaces on both edges
      splitstr_tovec <- function(Net_str_i) stringr::str_trim(unlist(strsplit(Net_str_i, sep, fixed=TRUE)), side = "both")
      # Turn a vector of character IDs into integer row numbers
      getRowsfromIDs <- function(NetIDvec) as.integer(sapply(NetIDvec, function(x) which(IDs %in% x)))
      # Turn any vector of IDs into a vector of length Kmax, filling remainder with trailing NA's
      makeKmaxIDvec <- function(NetIDVec) c(as.integer(NetIDVec), rep_len(NA_integer_, Kmax - length(NetIDVec)))

      Net_str <- as.character(data[,NETIDnode])
      NetIDs_l <- lapply(Net_str, splitstr_tovec) # Get list of n NET ID (character) vectors from NETIDnode
      NetRows_l <- NetIDs_l

      # if !is.null(IDnode), get the network row #s from IDs:
      if (!is.null(IDnode)) {
        IDs <- as.vector(data[, IDnode])
        NetRows_l <- lapply(NetIDs_l, getRowsfromIDs)
      }
      # NetVec_l <<- NetRows_l (NOT USED)
      # Make an array (n x Kmax) of network rows (filling remainder of each row with NA's)
      NetInd_k <- t(sapply(NetRows_l, makeKmaxIDvec))
      # NetInd_k <<- t(sapply(NetRows_l, makeKmaxIDvec))
      nF <- as.integer(.rowSums(! is.na(NetInd_k), m = nobs, n = Kmax))
      # nF <<- as.integer(.rowSums(! is.na(NetInd_k), m = nobs, n = Kmax))
      return(list(nF = nF, NetInd_k = NetInd_k)) # invisible(self)
    },

    mat.nF = function(nFnode) {
      assert_that(is.string(nFnode))
      nF <- as.matrix(nF)
      colnames(nF) <- nFnode
      # colnames(nF) <- "nF"
      nF
    }
  ),
  active = list(
    placeholder = function() {}
    # getNetInd = function() NetInd_k
  )
)

## ---------------------------------------------------------------------
# DETECTING VECTOR TYPES
## ---------------------------------------------------------------------
# sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
is.integerish <- function (x) {
  is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))
}
detect.col.types <- function(sVar_mat){
  sVartypes <- gvars$sVartypes
  as.list(apply(sVar_mat, 2,
                function(vec) {
                  vec_nomiss <- vec[!gvars$misfun(vec)]
                  nvals <- length(unique(vec_nomiss))
                  if (nvals <= 2L) {
                    sVartypes$bin
                  } else if ((nvals <= gvars$maxncats) && (is.integerish(vec_nomiss))) {
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
normalize_matsVar <- function(sVar_mat) {
  apply(sVar_mat, 2, normalize_sVar)
}
# Define bin cutt-offs for continuous x:
define.intervals <- function(x, nbins = gvars$nbins, bin_bymass = TRUE) {
  x <- x[!gvars$misfun(x)]  # remove missing vals
  nvals <- length(unique(x))
  if (nvals < nbins) nbins <- nvals # if nbins is too high, for ordinal, set nbins to n unique obs

  if (abs(max(x) - min(x)) > gvars$tolerr) {  # when x is not constant
    intvec <- seq.int(from = min(x), to = max(x), length.out = (nbins + 1)) # interval type 1: bin x by equal length intervals of 0-1
  } else {  # when x is constant, force the smallest possible interval to be at least [0,1]
    intvec <- seq.int(from = min(0L, min(x)), to = max(1L, max(x)), length.out = (nbins + 1))
  }
  if (bin_bymass) intvec <- quantile(x = x, probs = normalize(intvec)) # interval type 2: bin x by mass (quantiles of 0-1 intvec as probs)
  return(intvec) # return(list(intbylen = intvec, intbymass = intvecq))
}
# Turn any x into ordinal (1, 2, 3, ..., nbins) for a given interval cutoffs (length(intervals)=nbins+1)
make.ordinal <- function(x, intervals) findInterval(x = x, vec = intervals, rightmost.closed = TRUE)

# Remove the column naming for dummies_mat or keep for reference?
# Make dummy indicators for continuous x (sA[j])
make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms) {
# make.bins_mtx_1 = function(x.ordinal, nbins, bin.nms) { # Make dummy indicators for continuous x (sA[j])
  # Approach 1: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
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
# make.bins_mtx_2A = function(x.ordinal, self) { # (4) Approach 2A using model.matrix (a bit slower than manual)
#   # Creates B_j that jumps to 1 and then back to 0 (degenerate) and uses category 1 as ref
#   # First need to make a factor from xcat; adds a column of 1's (ref)
#   n <- length(x.ordinal)
#   cats = self$cats
#   x.ordinalf <- factor(x.ordinal);
#   dummies_mat <- model.matrix(~x.ordinalf);
#   dummies_mat
# }
# make.bins_mtx_2B = function(x.ordinal, self) { # (4) Approach 2B (manual vs. of 2A, still slower than 1)
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
#' \item{nodes} ...
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
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.flag
##' @export
DatNet <- R6Class(classname = "DatNet",
  portable = TRUE,
  class = TRUE,
  public = list(
    Kmax = integer(),          # max n of Friends in the network
    nodes = list(),            # names of the nodes in the data (Anode, Wnodes, nFnode, etc..)
    netind_cl = NULL,          # class NetIndClass object holding $NetInd_k network matrix
    VarNodes = NULL,           # Variable names (Wnodes or Anode)
    addnFnode = FALSE,         # Flag to add Fnode to output mat / df
    # misValRepl = FALSE,      # Replace missing network VarNode values (due to |nF| < Kmax) with gvar$misXreplace?
    # usemisval = NA,          # (NO LONGER FUNCTIONAL) Replacement value used for missing net covars (default is NA)
    ord.sVar = NULL,           # Ordinal (cat) transform for continous sVar

    active.bin.sVar = NULL,    # name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    # dat.bin.sVar = NULL,     # (MOVED TO AN ACT BIND) points to self$mat.bin.sVar
    mat.bin.sVar = NULL,       # temp storage mat for bin indicators on currently binarized continous sVar (from self$active.bin.sVar)

    mat.netVar = NULL,         # NOT DONE mat of network VarNode vals (+ VarNode itself) for each node in VarNodes
    # dat.netVar = NULL,       # (MOVED TO ACT BIND) df of network node values (+ node column itself) for all nodes in VarNodes

    Odata = NULL,              # data.frame used for creating the summary measures in mat.sVar, saved each time make.sVar called
    mat.sVar = NULL,           # Matrix storing all evaluated sVars, with named columns
    # dat.sVar = NULL,         # (MOVED TO ACT BIND) Matrix of summary meaure values for each summary measure expression in sVar

    sVar.object = NULL,        # Define_sVar object that contains / evaluates sVar expressions    
    
    type.sVar = NULL,          # named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    norm.c.sVars = FALSE,      # flag = TRUE if want to normalize continous covariates
    # nbins = list(),          # (MOVED TO ACT BIND) named list of total nbins for each "contin" sVar[i], can be overridden
    cbin_intrvls = list(),     # named list of bin cutoffs for each "contin" sVar[i], can be overridden

    # #todo 37 (DatNet) +0: NEED TO SORT OUT WHEN nOdata is needed....
    # nobs = NA_integer_,      # n of samples in the OBSERVED (original) data
    nOdata = NA_integer_,      # n of samples in the OBSERVED (original) data


    # replace with new version (VarNodes are no longer needed)
    # initialize = function(netind_cl, nodes, addnFnode = FALSE, ...) {
    initialize = function(netind_cl, nodes, VarNodes, addnFnode = FALSE, ...) {

      assert_that(!is.null(VarNodes)); assert_that(is.list(VarNodes)||is.character(VarNodes))
      self$VarNodes <- VarNodes 
      assert_that(is.flag(addnFnode))
      
      self$netind_cl <- netind_cl
      self$Kmax <- netind_cl$Kmax
      self$nodes <- nodes
      self$addnFnode <- addnFnode
      
      # NO LONGER FUNCTIONAL, AS misValRepl is now defined inside summary measure constructors:
      # assert_that(is.flag(misValRepl))
      # self$usemisval <- ifelse(misValRepl, gvars$misXreplace, gvars$misval)
      # NOTE: Eliminated this step in favor of cmisval being assigned the target value. This might be changed in the future.
      # Replace missing vals with gvars$misXreplace:
      # self$misValRepl <- misValRepl
      # if (self$misValRepl) self$fixmiss_netVar()

      invisible(self)
    },

    # New way of defining sA means no more need to always construct entire  mat.netVar
    # make.netVar.depr = function(Odata) {
    #   assert_that(is.data.frame(Odata))
    #   self$nOdata <- nrow(Odata)

    #   names.netVar <- netvar(unlist(self$VarNodes), (0L:self$Kmax))
    #   if (self$addnFnode) names.netVar <- c(names.netVar, self$nodes$nFnode)

    #   # Build network vectors: (W, W_netF_1, ..., W_netF_k) for each W in Wnodes
    #   self$mat.netVar <- matrix(0L, nrow = self$nOdata, ncol = length(names.netVar))
    #   colsperVar <- (self$Kmax + 1)
    #   for (idxVar in seq(self$VarNodes)) {
    #     Varnode <- self$VarNodes[[idxVar]]

    #     # Choose value to use for missing network covar: NA or 0
    #     # cmisval <- ifelse(self$misValRepl, gvars$misXreplace, gvars$misval)

    #     # THIS CAN BE REPLACED BY Define_sVar$parse with expr Varnode[[0:Kmax]]:
    #     self$mat.netVar[, ((idxVar - 1) * colsperVar + 1) : (idxVar * colsperVar)] <-
    #       .f.allCovars(k = self$Kmax, NetInd_k = self$netind_cl$NetInd_k, Var = Odata[,Varnode], VarNm = Varnode, misval = self$usemisval)
    #   }
    #   if (self$addnFnode) self$mat.netVar[, length(names.netVar)] <- as.matrix(Odata[, self$nodes$nFnode, drop=FALSE])
    #   colnames(self$mat.netVar) <- names.netVar
    # },

    # **********************
    # Define summary measures sVar
    # **********************
    # type.sVar acts as a flag: only detect types when !is.n, addnFnode = TRUEu
    make.sVar = function(Odata, sVar.object = NULL, type.sVar = NULL, norm.c.sVars = FALSE) {
      assert_that(is.data.frame(Odata))
      self$nOdata <- nrow(Odata)
      self$Odata <- Odata

      if (is.null(sVar.object)) {
        message("Not Implemented. To Be replaced with netVar construction when sVar.object is null...")
        stop()
      }

      self$sVar.object <- sVar.object
      if (self$addnFnode) { nFnode <- self$nodes$nFnode } else { nFnode <- NULL }

      self$mat.sVar <- sVar.object$get.mat.sVar(data.df = Odata, netind_cl = self$netind_cl, addnFnode = nFnode)

      # self$mat.sVar <- sVar.object$get.mat.sVar(data.df = Odata, netind_cl = self$netind_cl,
      #                                           misXreplace = self$usemisval, addnFnode = nFnode)

      # below was replaced with new sVar names that aren't nec. part of netVar:
      # assert_that(all(names.sVar %in% self$names.netVar))
      # n.sVar <- length(names.sVar)
      # mat.netVar <- self$mat.netVar
      # self$mat.sVar <- self$mat.netVar[, names.sVar]

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

    # *** TO DO ***: Make sure that a categorical var is only binned when ncats > nbins
    # *** TO DO ***: Make sure cat sVar doesn't get degenerate intervals (bins) (repeating same cut-offs)
    # Add arg (overwrite = FALSE)
    # if (!overwrite) and length(self$cbin_intrvls) > 0 => return(self$cbin_intrvls), otherwise recalculate, overwrite self$cbin_intrvls and return(self$cbin_intrvls)
    # Define the bin cut-off intervals for continous sVars
    # cbin_intrvls arg acts as a flag: only detect intervals when the arg is missing
    # otherwise can pass list(sVar = NA, ...) or as a default vector: c(0, 0.1, ..., 1)
    def_cbin_intrvls = function(cbin_intrvls, nbins = gvars$nbins, bin_bymass = TRUE) {
      # 1) select cont. sVars, proceed if length(c.sVars) > 0
      names.c.sVar <- self$names.c.sVar
      if (length(names.c.sVar) == 0L) {
        self$cbin_intrvls <- list()
        return(invisible(self))
      }
      # 2) for each c.sVars define (nbins+1) bin cutoff points
      # Unless its ordinal and ncats < nbins then define nbins = ncats
      if (missing(cbin_intrvls)) {
        cbin_intrvls <- vector(mode = "list", length = length(names.c.sVar))
        names(cbin_intrvls) <- names.c.sVar
        for (idx in seq_along(names.c.sVar)) {
          cbin_intrvls[[idx]] <- self$detect.sVar.intrvls(name.sVar = names.c.sVar[idx], nbins = nbins, bin_bymass = bin_bymass)
        }
        self$cbin_intrvls <- cbin_intrvls
      }
      # **** THE ROBUSTNESS CHECKS ARE NOT FINISHED **** 
      # 3) if !is.null(cbin_intrvls), check that cbin_intrvls is either a list of length(c.sVars) and that names match or 
      # 4) that cbin_intrvls is a vector & is.numeric(cbin_intrvls) & 
      # (max(cbin_intrvls) >= 1 && min(cbin_intrvls) <= 0)
      # or: (max(cbin_intrvls) >= max(sVar) && min(cbin_intrvls) <= min(norm(sVar)))
      if (!missing(cbin_intrvls)) {  # allows overriding of self$cbin_intrvls
        # ... check that names match before overwriting
        if (is.list(cbin_intrvls)) {
          self$cbin_intrvls <- cbin_intrvls
        } else if (is.numeric(cbin_intrvls)){
          self$cbin_intrvls <- rep(list(cbin_intrvls), length(names.c.sVar))
          names(self$cbin_intrvls) <- names.c.sVar
        }
      }
      invisible(self)
    },

    # #todo 18 (DatNet, DatNet.sWsA) +0: (OPTIONAL) ENABLE ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO DatNet
    # #todo 25 (DatNet, add_det) +0: Need to save det node flags as a separate mat, can't add them to sVar since all sVars will be automatically added to A ~ predictors
    add_deterministic = function(Odata, userDETcol) {
      determ.g_user <- as.vector(Odata[,userDETcol]) # get deterministic As for the entire network of each unit (set by user)
      # determ.gvals_user <- Odata[,AnodeDET] # add values to be assigned to deterministic nodes (already have: netA)
      determ_cols_user <- .f.allCovars(k = self$Kmax, NetInd_k = self$netind_cl$NetInd_k, Var = determ.g_user,
                                        VarNm = "determ.g_true", misval = gvars$misval)
      # determ_cols <- (determ_cols_user | determ_cols_Friend)
      determ_cols <- determ_cols_user
      # THIS IS A VERY BAD WAY TO DO THAT, REMOVING
      # self$mat.sVar <- cbind(determ_cols, self$mat.sVar)
      determ_cols_type <- as.list(rep_len(gvars$sVartypes$bin, ncol(determ_cols)))
      names(determ_cols_type) <- colnames(determ_cols)
      self$type.sVar <- c(determ_cols_type, self$type.sVar)
      invisible(self)
    },
    # Replaces all misval values in mat.netVar / mat.sVar with gvars$misXreplace
    fixmiss_netVar = function() {
      self$mat.netVar[gvars$misfun(self$mat.netVar)] <- gvars$misXreplace
      invisible(self)
    },
    fixmiss_sVar = function() {
      self$mat.sVar[gvars$misfun(self$mat.sVar)] <- gvars$misXreplace
      invisible(self)
    },

    # --------------------------------------------------
    # Methods for directly handling one continous/categorical sVar in self$mat.sVar;
    # No checking of incorrect input is performed, use at your own risk!
    # MIGHT NEED TO MOVE TOP METHODS (discretize.sVar & binirize.sVar) TO DatNet.sWsA (see note in h.fit())
    # --------------------------------------------------
    norm.sVar = function(name.sVar) { normalize_sVar(self$dat.sVar[, name.sVar]) },  # return normalized 0-1 sVar
    set.sVar = function(name.sVar, new.sVar) { self$mat.sVar[, name.sVar] <- new.sVar },
    get.sVar = function(name.sVar) { self$dat.sVar[, name.sVar] },
    set.sVar.type = function(name.sVar, new.type) { self$type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { self$type.sVar } else { self$type.sVar[[name.sVar]] } },
    
    # *** NOTE *** When sVar is cat might be better to set bin_bymass = FALSE to avoid collapsing of categories for sVar
    # detect.sVar.intrvls = function(name.sVar, nbins = gvars$nbins, bin_bymass = FALSE) {
    detect.sVar.intrvls = function(name.sVar, nbins = gvars$nbins, bin_bymass = TRUE) {
      int <- define.intervals(x = self$dat.sVar[, name.sVar], nbins = nbins, bin_bymass = bin_bymass)
      if (length(unique(int)) < length(int)) {
        message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+% (length(int)-1) %+% " to " %+% (length(unique(int))-1) %+% " due to too few obs.")
        print("old intervals: "); print(int)
        int <- unique(int)
        print("new intervals: "); print(int)
      }
      return(int)
    },

    set.sVar.intrvls = function(name.sVar, new.intrvls) { self$cbin_intrvls[[name.sVar]] <- new.intrvls },
    get.sVar.intrvls = function(name.sVar) { if (missing(name.sVar)) { self$cbin_intrvls } else { self$cbin_intrvls[[name.sVar]] } },
    # Need to find a way to over-ride nbins for categorical vars (allowing it to be set to more than gvars$maxncats)!
    nbins.sVar = function(name.sVar) { if (missing(name.sVar)) { self$all.nbins } else { length(self$get.sVar.intrvls(name.sVar)) - 1 } },
    # Return names of bin indicators for sVar:
    bin.nms.sVar = function(name.sVar) { name.sVar%+%"_"%+%"B."%+%(1L:self$nbins.sVar(name.sVar)) }, 
    # create a vector of ordinal (categorical) vars out of cont. sVar vector:
    discretize.sVar = function(name.sVar, intervals = self$get.sVar.intrvls(name.sVar)) {  
      self$ord.sVar <- make.ordinal(x = self$dat.sVar[, name.sVar], intervals = intervals)
      invisible(self$ord.sVar)
    },
    # return matrix of bin indicators for ordinal
    binirize.sVar = function(name.sVar, intervals = self$get.sVar.intrvls(name.sVar), nbins = self$nbins.sVar(name.sVar), bin.nms = self$bin.nms.sVar(name.sVar)) {
      self$mat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$discretize.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms)
      self$active.bin.sVar <- name.sVar
      invisible(self$mat.bin.sVar)
    }
  ),

  active = list(
    names.netVar = function() { colnames(self$dat.netVar) },
    names.sVar = function() { colnames(self$dat.sVar) },
    names.c.sVar = function() { names(self$type.sVar[self$type.sVar %in% gvars$sVartypes$cont]) }, # names of cont sVars
    all.nbins = function() { lapply(self$cbin_intrvls, function(int) length(int)-1) },
    ncols.netVar = function() { length(self$names.netVar) },
    ncols.sVar = function() { length(self$names.sVar) },

    dat.netVar = function() { self$mat.netVar },
    dat.sVar = function() { self$mat.sVar },
    dat.bin.sVar = function() { self$mat.bin.sVar },

    emptydat.netVar = function() {self$mat.netVar <- NULL },   # wipe out mat.netVar
    emptydat.sVar = function() { self$mat.sVar <- NULL }       # wipe out mat.sVar
  ),
  private = list(
    # mat.netVar = NULL,
    # mat.sVar = NULL,
    # mat.bin.sVar = NULL,
    placeholder = list()
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
#' \item{datnetW} ...
#' \item{datnetA} ...
#' \item{YnodeVals} ...
#' \item{det.Y} ...
#' \item{Kmax} ...
#' }
#' Methods for combining, subsetting, discretizing & binirizing of summary measures in sW & sA. Inherits from DatNet class. 
#' The combined dataset of all (sW, sA) summary measures returned by DatNet is stored as a matrix under self$mat.sVar
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.flag
##' @export
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
      print("entered DatNet.sWsA constructor")
      assert_that("DatNet" %in% class(datnetW))
      assert_that("DatNet" %in% class(datnetA))

      self$datnetW <- datnetW
      self$datnetA <- datnetA
      self$Kmax <- datnetW$Kmax
      self$nodes <- datnetW$nodes
      self$netind_cl <- datnetW$netind_cl

      if (!missing(YnodeVals)) {
        if (missing(det.Y)) det.Y <- rep_len(FALSE, length(YnodeVals))
        self$noNA.Ynodevals <- YnodeVals  # Adding actual observed Y as protected (without NAs)
        self$YnodeVals <- YnodeVals
        self$YnodeVals[det.Y] <- NA       # Adding public YnodeVals & setting det.Y values to NA
        self$det.Y <- det.Y
      }

      invisible(self)
    },

    # x) Note that if this feature is to be open to the user, eval() has to be done in the parent.frame of the calling function, not baseenv().
    # x) Same with summary measures: need to eval them in the calling environment (in addition to the envir of data.frame(netW,netA))
    evalsubst = function(subsetexpr) { # Eval the expression (in the environment of the data.frame "data" + global constants "gvars"):      
      if (is.logical(subsetexpr)) {
        return(subsetexpr)
      } else {
        # REPLACING WITH env that is made of data.frames instead of matrices
        # # todo 31 (evalsubst) +0: NEED TO RE-WRITE ALL SUBSET EVALUATIONs SO THAT IT WORKS WITH MATRICES (can't use expressions for subset anymore)
        # #todo 36 (evalsubst) +0: Possible ideas, instead of sVar.name inside the expression !misfun(sum_1mAW2_nets), use !misfun(dat.sWsA[,sVar.name]) / !misfun(dat.sVar[,sVar.name]) or !misfun(dat.bin.sVar[,sVar.name])
        # Could do evaluation in a special env with a custom subsetting fun '[' that will dynamically find the write dataset that contains sVar.name (dat.sVar or dat.bin.sVar) and will return sVar vector
        eval.env <- c(data.frame(self$dat.sWsA), data.frame(self$dat.bin.sVar), as.list(gvars))
        # eval.env <- c(self$dat.sWsA, self$dat.bin.sVar, as.list(gvars))
        res <- try(eval(subsetexpr, envir = eval.env, enclos = baseenv())) # to evaluate vars not found in data in baseenv()      
        return(res)
        # Older versions:
        # res <- try(eval(subsetexpr, envir = c(self$datnetW$dat.sVar, self$datnetA$dat.sVar, self$datnetA$dat.bin.sVar, as.list(gvars)), enclos = baseenv())) # to evaluate vars not found in data in baseenv()
        # res <- try(eval(subst_call, envir = data, enclos = parent.frame())) # to evaluate vars not found in data in parent.frame()
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
      # LATEST VERSION (data.frame sWsA is saved as a self$dat.sWsA):
      dfsel <- self$dat.sWsA[rowsubset, sel.sWsA, drop = FALSE]
      if (length(sel.binsA)>0) {
        dfsel <- cbind(dfsel, dat.bin.sVar[rowsubset, sel.binsA, drop = FALSE])
      }

      # OLDER VERSIONs. This version is probably more memory eff, but takes longer and is very prone to bugs
      # Allows datnetW to to be stored with only nrow = nobs(Odata), while datnetA has nrow = p*nobs(Odata) under g_star
      # datnetW <- self$datnetW
      # datnetA <- self$datnetA
      # sel.sW <- TRUE
      # sel.sA <- TRUE
      # dat.bin.sVar <- self$datnetA$dat.bin.sVar
      # selrow1 <- selrow2 <- selrow3 <- rowsubset
      # if (length(rowsubset) > nrow(dat.bin.sVar))     selrow1 <- rowsubset[1:nrow(dat.bin.sVar)]
      # if (length(rowsubset) > nrow(datnetA$dat.sVar)) selrow2 <- rowsubset[1:nrow(datnetA$dat.sVar)]
      # if (length(rowsubset) > nrow(datnetW$dat.sVar)) selrow3 <- rowsubset[1:nrow(datnetW$dat.sVar)]
      # dfsel <- cbind(dat.bin.sVar[selrow1, sel.binsA, drop = FALSE],
      #                 datnetA$dat.sVar[selrow2, sel.sA, drop = FALSE],
      #                 datnetW$dat.sVar[selrow3, sel.sW, drop = FALSE])
      # VERSION 2: DOESN'T WORK when dat.bin.sVar is empty (data.frame()) -> can't cbind 0 rows with N>0 rows. VERY BUGGY
      # Select the columns, cbind together, then select the rows once:
      # dfsel <- cbind(dat.bin.sVar[, sel.binsA, drop = FALSE],
      #                 datnetA$dat.sVar[, sel.sA, drop = FALSE],
      #                 datnetW$dat.sVar[, sel.sW, drop = FALSE])
      # dfsel <- dfsel[rowsubset, , drop = FALSE]
      return(dfsel)
    },

    get.outvar = function(rowsubset = TRUE, var) {
      if (var %in% self$names.sWsA) {
        self$dat.sWsA[rowsubset, var]
      } else if (var %in% colnames(self$dat.bin.sVar)) {
        self$dat.bin.sVar[rowsubset, var]
      } else if ((var %in% self$nodes$Ynode) && !is.null(self$YnodeVals)) {
        self$YnodeVals[rowsubset]
      } else {
        stop("outcome variable not found")
      }
    },

    copy.cbin.intrvls = function() {
      self$type.sVar <- c(self$datnetW$type.sVar, self$datnetA$type.sVar)
      self$cbin_intrvls <- c(self$datnetW$cbin_intrvls, self$datnetA$cbin_intrvls)
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
      # Copy variable detected types and bin interval definitions from the observed data classes (datnetW, datnetA) to self:
      self$copy.cbin.intrvls()

      if (is.null(f.g_name)) {  # set df.sWsA to observed data (sW,sA) if g.fun is.null
        df.sWsA <- cbind(datnetW$dat.sVar, datnetA$dat.sVar) # assigning summary measures as data.frames:
      } else {  # need to sample A under f.g_name (gstar or known g0), possibly re-evaluate sW from O.datnetW

        Odata <- datnetW$Odata

        # With new structure there is only one sW, which is built only once in datnetW
        # if (!is.null(sW.object)) { # re-creating sW, e.g., when using sW.gstar different from sW.g0
        #   self$datnetW <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes, VarNodes = self$nodes$Wnodes, addnFnode = TRUE)
        #   # self$datnetW <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes, VarNodes = self$nodes$Wnodes, addnFnode = TRUE, misValRepl = TRUE)
        #   self$datnetW$make.sVar(Odata = Odata, sVar.object = sW.object)
        #   self$datnetW$fixmiss_sVar()
        # }

        # Should we save new datnetA.gstar as self$datnetA (removing a pointer to O.datnetA)? Not for now.
        datnetA.gstar <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes, VarNodes = self$nodes$Anode)
        # datnetA.gstar <- DatNet$new(netind_cl = datnetW$netind_cl, nodes = self$nodes, VarNodes = self$nodes$Anode, misValRepl = TRUE)

        df.sWsA <- matrix(nrow = (nobs * p), ncol = (datnetW$ncols.sVar + datnetA$ncols.sVar))  # pre-allocate result matx sWsA
        colnames(df.sWsA) <- self$names.sWsA

        # df.sWsA <- data.frame(df.sWsA) # NOT SURE WE STILL WANT TO CONVERT THIS TO df?

        for (i in seq_len(p)) {

          # *** f.g_name can only depend on covariates in datnetW$dat.sVar ***
          A.gstar <- f.gen.A.star(self$Kmax, datnetW$dat.sVar, f.g_name, f.g_args)
          # Avec.df <- data.frame(Anode = f.gen.A.star(self$Kmax, datnetW$dat.sVar, f.g_name, f.g_args), stringsAsFactors = FALSE)
          # colnames(Avec.df)[1] <- self$nodes$Anode
          Odata[, self$nodes$Anode] <- A.gstar # replace A under g0 in Odata with A^* under g.star:
          datnetA.gstar$make.sVar(Odata = Odata, sVar.object = sA.object) # create new summary measures sA (under g.star)
          # Assiging summary measures to one output data:
          df.sWsA[((i - 1) * nobs + 1):(nobs * i), ] <- cbind(datnetW$dat.sVar, datnetA.gstar$dat.sVar)[,]
        }
      }

      self$mat.sVar <- df.sWsA

      invisible(self)
    }
  ),
  active = list(
    dat.sWsA = function() { self$mat.sVar }, # NO LONGER NEEDED, REMOVE
    # df.sWsA = function() { self$mat.sVar }, # REMOVED / was used for compatibility with older functions
    names.sWsA = function() { c(self$datnetW$names.sVar, self$datnetA$names.sVar) },
    nobs = function() { nrow(self$dat.sWsA) },
    noNA.Ynodevals = function(noNA.Yvals) {
      if (missing(noNA.Yvals)) return(private$protected.YnodeVals)
      else private$protected.YnodeVals <- noNA.Yvals
    }
  ),
  private = list(
    protected.YnodeVals = NULL, # Actual observed values of the binary outcome (Ynode), along with deterministic vals
    placeholder = list()
  )
)