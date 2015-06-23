#-----------------------------------------------------------------------------
# MAJOR CHANGES TO DatNet.sWsA CLASS STRUCTURE (06/18/2015):
#-----------------------------------------------------------------------------
# * gen_sWsA_dat() renamed to make.df.sWsA() and moved inside DatNet.sWsA class;
# * The observed data is now always stored as $datnetW,$datnetA fields in every DatNet.sWsA object (under g_0 or g_star);
# * DatNet.sWsA$new(datnetW,datnetW) only saves obs data objects datnetW and datnetW as fields;
# * DatNet.sWsA$new() has to be called only twice: once for g_0 (g_N) and once for g_star;
# * $make.df.sWsA() must be called to create $df.sW.sA ($df.sVar) - a data.frame of COMBINED sWsA;
# * All binning / interval methods MUST BE called on DatNet.sWsA (NOT $datnetA, $datnetW) (inherited from DatNet);

# * Both copies of DatNet.sWsA are storing datnetA/datnetW by reference - same copy;
# * Changing datnetA/datnetW in one copy of DatNet.sWsA will result them being changed in the other copy of DatNet.sWsA as well.

# * MOST IMPORTANTLY we can now get rid of the MC sim loop for evaling psi_n. Just use already sampled DatNet.gstar dataset and evaluate psi_n only once.



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

## ---------------------------------------------------------------------
# Class for managing and generating network/summary matrices netVar/sVar where Var is Wnodes or Anode
## ---------------------------------------------------------------------
# TO DO: 
# x) Finish sVar generation (from expressions)
# x) See if some of new() args can be dropped: (Odata, NetInd_k, Kmax, nodes, VarNodes, AddnFnode = FALSE, misValRepl = FALSE, ...) {
# x) See if new() can be made to sample Odata[,Varnode] if f.g_name is supplied... 
  # This would allow making DatNet fully closed in, that is, its the only class where netVar and sVar data are contructed, given Odata and f.g_name...

#-----------------------------------------------------------------------------
# - Subsetting / Delayed eval
#-----------------------------------------------------------------------------
  # x) Moved delayed eval. to be inside DatNet.sWsA$evalsubst
  # x) Delayed evaluation of subsets implies that data / newdata have to be data.frames at the time of subseting. 
  # x) See implementation of subset.data.frame and subset.matrix - there might be a way around it even for data=as.matrix(data)

#-----------------------------------------------------------------------------
# TO DO: 
# - Replace all node name references (columns) with indices? See ?subset.data.frame:
  # nl <- as.list(seq_along(data.df))
  # names(nl) <- names(data.df)
  # eval(substitute(node), nl, parent.frame()) -> replace parent.frame() with a reasonable calling envir
## ---------------------------------------------------------------------



## ---------------------------------------------------------------------
# NOT FINISHED
# Very simple class to hold and create NetInd_k, the matrix of network connection indices in observed data, dim = (nobs x Kmax)
## ---------------------------------------------------------------------
NetIndClass <- R6Class("NetIndClass",
  class = FALSE,
  portable = FALSE,
  public = list(
    NetInd_k = matrix(),
    nobs = NA_integer_,
    Kmax = NA_integer_,
    initialize = function(Odata, Kmax) {
      nobs <<- nrow(Odata)
      Kmax <<- Kmax
      NetInd_k <<- matrix(0L, nrow = nobs, ncol = Kmax) # pre-allocate NetInd_k mat
      # makeNetInd(Odata = Odata)
    },
    assignNetInd = function(NetInd_k) {
      NetInd_k[] <<- NetInd_k
      # or below?
      # NetInd_k[,] <<- NetInd_k
    },
    makeNetInd = function(Odata) {
      # method for constructing NetInd_k
      # NetInd_k <<- matrix(0L, nrow = nobs, ncol = Kmax)
      # ....
      # ....
    }
  ),
  active = list(
    getNetInd = function() NetInd_k
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
normalize = function(x) {
  if (abs(max(x) - min(x)) > gvars$tolerr) { # Normalize to 0-1 only when x is not constant
    return((x - min(x)) / (max(x) - min(x)))
  } else {  # What is the thing to do when x constant? Set to abs(x), abs(x)/x or 0???
    return(x)
  }
}
normalize_sVar = function(sVar_vec) {
  nonmiss_idx <- !gvars$misfun(sVar_vec)
  if (sum(nonmiss_idx) > 0) {
    sVar_vec[nonmiss_idx] <- normalize(sVar_vec[nonmiss_idx])
  }
  sVar_vec
}
normalize_matsVar = function(sVar_mat) {
  apply(sVar_mat, 2, normalize_sVar)
}
define.intervals = function(x, nbins = gvars$nbins, bin_bymass = TRUE) {  # Define bin cutt-offs for continuous x.
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
make.ordinal = function(x, intervals) { 
  findInterval(x = x, vec = intervals, rightmost.closed = TRUE)
}
# Remove the column naming for dummies_mat or keep for reference?
make.bins_mtx_1 = function(x.ordinal, nbins, bin.nms) { # Make dummy indicators for continuous x (sA[j])
# make.bins_mtx_1 = function(x.ordinal, nbins, bin.nms) { # Make dummy indicators for continuous x (sA[j])
  # Approach 1: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
  n <- length(x.ordinal)
  cats <- (1:nbins)
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
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.flag
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
    ord.sVar = NULL,           # ordinal transform for regression outcome sVar

    # bin_mat.sVar = NULL,     # (ACTIVE BINDING) matrix of indicators for binarized regression outcome sVar
    df.bin.sVar = NULL,		   # data.frame of indicators for binarized regression outcome sVar

    dat.netVar = NULL,         # matrix of network node values (+ node column itself) for all nodes in VarNodes
    # df.netVar = NULL,        # (ACTIVE BINIDNG) df of network node values (+ node column itself) for all nodes in VarNodes

    # dat.sVar = NULL,         # (ACTIVE BINIDNG) Matrix of summary meaure values for each summary measure expression in sVar
    df.sVar = NULL,            # data.frame of summary meaure values for each summary measure expression in sVar

    type.sVar = NULL,          # named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    norm.c.sVars = FALSE,      # flag = TRUE if want to normalize continous covariates
    # nbins = list(),          # (ACTIVE BINDING) named list of total nbins for each "contin" sVar[i], can be overridden
    cbin_intrvls = list(),     # named list of bin cutoffs for each "contin" sVar[i], can be overridden
    
    # NEED TO SORT OUT WHEN nOdata is needed....
    # nobs = NA_integer_,      # n of samples in the OBSERVED (original) data
    nOdata = NA_integer_,      # n of samples in the OBSERVED (original) data

    initialize = function(Odata, NetInd_k, Kmax, nodes, VarNodes, AddnFnode = FALSE, misValRepl = FALSE, ...) {
      assert_that(!is.null(VarNodes)); assert_that(is.list(VarNodes)||is.character(VarNodes))
      self$VarNodes <- VarNodes # self$VarNodes <- nodes$Wnodes # or: # self$VarNodes <- nodes$Anode
      # a better alternative is: # self$VarNodes <- nodes[[VarNodeType]] # VarNodeType <- "Wnodes"; VarNodeType <- "Anode"
      assert_that(is.count(Kmax))
      assert_that(is.flag(AddnFnode))
      assert_that(is.flag(misValRepl))
      self$NetInd_k <- NetInd_k   # include or not?
      self$Kmax <- Kmax
      self$nodes <- nodes
      assert_that(is.data.frame(Odata)) 
      self$nOdata <- nrow(Odata)
      self$AddnFnode <- AddnFnode
      self$misValRepl <- misValRepl
      names.netVar <- netvar2(unlist(self$VarNodes), (0L:self$Kmax))
      if (self$AddnFnode) names.netVar <- c(names.netVar, self$nodes$nFnode)
      cmisval <- ifelse(self$misValRepl, gvar$misXreplace, gvars$misval) # Use gvars$misXreplace instead of the missing values. This might be changed in the future.
      # Build network vectors: (W, W_netF_1, ..., W_netF_k) for each W in Wnodes by PRE-ALLOCATING netW_full:
      self$dat.netVar <- matrix(0L, nrow = self$nOdata, ncol = length(names.netVar))
      colsperVar <- (self$Kmax + 1)
      for (idxVar in seq(self$VarNodes)) {
        Varnode <- self$VarNodes[[idxVar]]
        self$dat.netVar[, ((idxVar - 1) * colsperVar + 1) : (idxVar * colsperVar)] <- 
              .f.allCovars(k = self$Kmax, NetInd_k = self$NetInd_k, Var = Odata[,Varnode], VarNm = Varnode, misval = cmisval)
      }
      if (self$AddnFnode) self$dat.netVar[, length(names.netVar)] <- as.matrix(Odata[, self$nodes$nFnode, drop=FALSE])
      colnames(self$dat.netVar) <- names.netVar
      # Replace missing vals with gvars$misXreplace:
      # NOTE: Eliminated this step in favor of cmisval being assigned the target value. This might be changed in the future. 
      # if (self$misValRepl) self$fixmiss_netVar() 
      invisible(self)
    },
    # *** MAKE A PRIVATE METHOD ***
    # Define the type (class) of each summary measure: bin, cat or cont
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar)
    # otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    def_types_sVar = function(type.sVar = NULL) {
      if (is.null(type.sVar)) { # Detect the type of each sVar[i]: gvars$sVartypes$bin,  gvars$sVartypes$cat, gvars$sVartypes$cont
        self$type.sVar <- detect.col.types(self$df.sVar)
        
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
    # Normalize continuous sVars
    # This could be memory-costly. Alternative is a for loop over cols in self$df.sVar
    norm_c_sVars = function() {
      names.c.sVar <- self$names.c.sVar
      if (length(names.c.sVar) == 0L) return(invisible(self))

      if (self$norm.c.sVars && (length(names.c.sVar) > 0)) {
        for (name.c.sVar in names.c.sVar) {
          self$df.sVar[, name.c.sVar] <- normalize_sVar(self$df.sVar[, name.c.sVar])
        }
      }
      invisible(self)
    },
    # Define the summary measures (NOT IMPLEMENTED. )
    # For now can only specify summary measure as subsets of netVar columns
    # *** TO DO:
    # Apply the summary measure functions / expression to dat.netVar to OBTAIN df.sVar columns
    # Select only sW columns in hform_g0 and hfrom_gstar or use all?
    # type.sVar acts as a flag: only detect types when !is.null
    make_sVar = function(names.sVar = NULL, type.sVar = NULL, norm.c.sVars = FALSE) {
      if (is.null(names.sVar)) {
        names.sVar <- self$names.netVar
      }
      assert_that(is.character(names.sVar))
      # to be replaced with new sVar names that aren't nec. part of netVar:
      assert_that(all(names.sVar %in% self$names.netVar)) 
      n.sVar <- length(names.sVar)   

      dfnetVar <- self$df.netVar
      self$df.sVar <- dfnetVar[, names.sVar]
      # self$dat.sVar <- self$dat.netVar[, names.sVar]

      # MAKE def_types_sVar an active binding? calling self$def_types_sVar <- type.sVar assigns, calling self$def_types_sVar defines.
      self$def_types_sVar(type.sVar) # Define the type of each sVar[i]: bin, cat or cont
      # PUT IT OUTSIDE OF make_sVar as a separate public method?
      if (norm.c.sVars) { # normalize continuous and non-missing sVars, overwrite their columns in df.sVar with normalized [0,1] vals
        self$norm.c.sVars <- norm.c.sVars  
        self$norm_c_sVars()
      }
      invisible(self)
    },
    # MAKE INTO AN ACTIVE BINDING????  # YES. Add arg (overwrite = FALSE) 
    # if (!overwrite) and length(self$cbin_intrvls) > 0 => return(self$cbin_intrvls), otherwise recalculate, overwrite self$cbin_intrvls and return(self$cbin_intrvls)
    # *** TO DO ***: Make sure that a categorical var is only binned when ncats > nbins
    # Define the bin cut-off intervals for continous sVars
    # cbin_intrvls acts as a flag: only detect intervals when its missing
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
            cbin_intrvls[[idx]] <- define.intervals(x = self$df.sVar[, names.c.sVar[idx]], nbins = nbins, bin_bymass = bin_bymass)
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
    # (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sW:
    add_deterministic = function(Odata, userDETcol) {
      determ.g_user <- as.vector(Odata[,userDETcol]) # get deterministic As for the entire network of each unit (set by user)
      # determ.gvals_user <- Odata[,AnodeDET] # add values to be assigned to deterministic nodes (already have: netA)
      determ_cols_user <- .f.allCovars(k = self$Kmax, NetInd_k = self$NetInd_k, Var = determ.g_user, 
                                        VarNm = "determ.g_true", misval = gvars$misval)
      # determ_cols <- (determ_cols_user | determ_cols_Friend)
      determ_cols <- determ_cols_user
      self$df.sVar <- cbind(determ_cols, self$df.sVar)
      determ_cols_type <- as.list(rep_len(gvars$sVartypes$bin, ncol(determ_cols)))
      names(determ_cols_type) <- colnames(determ_cols)
      self$type.sVar <- c(determ_cols_type, self$type.sVar)
      invisible(self)
    },
    # REPLACE ALL misval values in sW with gvars$misXreplace (OR DO IT ON THE ACTUAL DATASET WHEN SUBSETTING IS PERFORMED)
    fixmiss_netVar = function() {
      self$dat.netVar[gvars$misfun(self$dat.netVar)] <- gvars$misXreplace
      invisible(self)
    },
    fixmiss_sVar = function() {
      self$df.sVar[gvars$misfun(self$df.sVar)] <- gvars$misXreplace
      invisible(self)
    },
    # --------------------------------------------------
    # Methods for directly handling one continous/categorical sVar in self$df.sVar;
    # No checking of incorrect input is performed, use at your own risk!
    # --------------------------------------------------
    norm.sVar = function(name.sVar) { normalize_sVar(self$df.sVar[, name.sVar]) },  # return normalized 0-1 sVar
    set.sVar = function(name.sVar, new.sVar) { self$df.sVar[, name.sVar] <- new.sVar },
    get.sVar = function(name.sVar) { self$df.sVar[, name.sVar] },
    set.sVar.type = function(name.sVar, new.type) { self$type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { self$type.sVar } else { self$type.sVar[[name.sVar]] } },
    detect.sVar.intrvls = function(name.sVar, nbins = gvars$nbins, bin_bymass = TRUE) { define.intervals(x = self$df.sVar[, name.sVar], nbins = nbins, bin_bymass = bin_bymass) },
    set.sVar.intrvls = function(name.sVar, new.intrvls) { self$cbin_intrvls[[name.sVar]] <- new.intrvls },
    get.sVar.intrvls = function(name.sVar) { if (missing(name.sVar)) { self$cbin_intrvls } else { self$cbin_intrvls[[name.sVar]] } },
    # Need to find a way to over-ride nbins for categorical vars (allowing to set to more than ncats)!
    nbins.sVar = function(name.sVar) { if (missing(name.sVar)) { self$all.nbins } else { length(self$get.sVar.intrvls(name.sVar)) - 1 } },
    bin.nms.sVar = function(name.sVar) { name.sVar%+%"_"%+%"B."%+%(1L:self$nbins.sVar(name.sVar)) }, # Return names of bins indicators for sVar
    make.ord.sVar = function(name.sVar, intervals = self$get.sVar.intrvls(name.sVar)) {  # create vector of ordinals out of cont. sVar
      self$ord.sVar <- make.ordinal(x = self$df.sVar[, name.sVar], intervals = intervals) 
      invisible(self$ord.sVar)
    },
    binirize.sVar = function(name.sVar, intervals = self$get.sVar.intrvls(name.sVar), nbins = self$nbins.sVar(name.sVar), bin.nms = self$bin.nms.sVar(name.sVar)) {  # return matrix of bin indicators for ordinal
      self$df.bin.sVar <- data.frame(make.bins_mtx_1(x.ordinal = self$make.ord.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms))
      invisible(self$df.bin.sVar)
    }
  ),
  active = list(
    names.netVar = function() { colnames(self$dat.netVar) },
    names.sVar = function() { colnames(self$df.sVar) },
    names.c.sVar = function() { names(self$type.sVar[self$type.sVar %in% gvars$sVartypes$cont]) }, # names of cont sVars
    all.nbins = function() { lapply(self$cbin_intrvls, function(int) length(int)-1) },
    ncols.netVar = function() { length(self$names.netVar) },
    ncols.sVar = function() { length(self$names.sVar) },
    df.netVar = function() { data.frame(self$dat.netVar) },    # convert matrix dat.netVar to a data.frame
    dat.sVar = function() { as.matrix(self$df.sVar) },		   # convert df df.sVar to a matrix
    # bin_mat.sVar = function() { as.matrix(self$df.bin.sVar) },
    emptydat.netVar = function() {self$dat.netVar <- NULL },   # wipe out dat.netVar
    emptydf.sVar = function() { self$df.sVar <- NULL }       # wipe out df.sVar
  ),
  private = list(
    # dat.netVar = NULL,
    # df.sVar = NULL,
    placeholder = list()
  )
)

## ---------------------------------------------------------------------
# Class for managing and generating all the summary datasets sWsA used in fitting \bar{h}^* or \bar{h}_0
# DatNet.sWsA is the only way to access data in the package 
# Gets passed on to SummariesModel functions: $fit(), $predict() and $predictAeqa()
## ---------------------------------------------------------------------
# TO DO: 
# x) Save data.frame(datbin_mat) (binarized continous sA) as a field and include it in eval_subset and subsetting functions
# x) MAKE DatNet.sWsA inherit from DatNet and DO NOT REPLICATE ANY FIELDS THAT ALREADY EXIST IN DatNet!
# x) df.sW.sA can be A LOT SMARTER and memory efficient. 
    # If we know the covariates that are needed it might do the subsetting $df.sW.sA(covars)
    # Then first select those covariate from datnetW$df.sVar and datnetA$df.sVar and return only what's needed.

DatNet.sWsA <- R6Class(classname = "DatNet.sWsA",
  inherit = DatNet,
  portable = TRUE,
  class = TRUE,
  public = list(
    datnetW = NULL,   # *** RENAME TO O.datnetW for clarity ***
    datnetA = NULL,   # *** RENAME TO O.datnetA for clarity ***
    # **********
    # df.sVar - (inherited): now stores the combined data.frame of cbind(data.frame(dat.sW), data.frame(dat.sA))
    # this way ALL the methods and active bindings of DatNet are valid in DatNet.sWsA for the combined data.frame
    # **********
    df.sW.sA = NULL, # assigned to self$df.sVar for compatibility with old functions
    # **********

    initialize = function(datnetW, datnetA, ...) {
      print("entered DatNet.sWsA constructor")
      assert_that("DatNet" %in% class(datnetW))
      assert_that("DatNet" %in% class(datnetA))
      self$datnetW <- datnetW
      self$datnetA <- datnetA
      invisible(self)
    },
    # x) Note that if this feature is to be open to the user, eval() has to be done in the parent.frame of the calling function, not baseenv().
    # x) Same with summary measures: need to eval them in the calling environment (in addition to the envir of data.frame(netW,netA))
    evalsubst = function(subsetexpr) { # Eval the expression (in the environment of the data.frame "data" + global constants "gvars"):      
      # print("DatNet.sWsA evalsubst....")
      # print("self$df.sW.sA"); print(head(self$df.sW.sA))
      # print("self$df.bin.sVar"); print(head(self$df.bin.sVar))
      res <- try(eval(subsetexpr, envir = c(self$df.sW.sA, self$df.bin.sVar, as.list(gvars)), enclos = baseenv())) # to evaluate vars not found in data in baseenv()      
      # Eval the expression in the environments of 3 data.frames + global constants "gvars"):
      # res <- try(eval(subsetexpr, envir = c(self$datnetW$df.sVar, self$datnetA$df.sVar, self$datnetA$df.bin.sVar, as.list(gvars)), enclos = baseenv())) # to evaluate vars not found in data in baseenv()
      # res <- try(eval(subst_call, envir = data, enclos = parent.frame())) # to evaluate vars not found in data in parent.frame()
      # old:# res <- try(eval(subst_call, envir = c(lapply(TD_vnames, I), node_func), enclos=anchor_evn))  # evaluate modified_call in the df namespace with custom '[' function  
     return(res)
    },
    get.df.sW.sA = function(rowsubset = TRUE, covars) { # return a data.frame with covars (predictors)
      df.bin.sVar <- self$df.bin.sVar
      sel.sWsA <- TRUE
      sel.binsA = NULL # columns to select from binned continuos var matrix (if it was previously constructed)
      if (!missing(covars)) {
        sel.sWsA <- colnames(self$df.sW.sA)[(colnames(self$df.sW.sA) %in% covars)]
        if (!is.null(df.bin.sVar)) {
          sel.binsA <- colnames(df.bin.sVar)[(colnames(df.bin.sVar) %in% covars)]
        }
      }
      # LATEST VERSION (data.frame sWsA is saved as a self$df.sW.sA):
      dfsel <- self$df.sW.sA[rowsubset, sel.sWsA, drop = FALSE]
      if (length(sel.binsA)>0) {
        dfsel <- cbind(df.bin.sVar[rowsubset, sel.binsA, drop = FALSE], dfsel)
      }
      # VERSION 1:
      # This version is probably more memory friendly, but takes longer time to run
      # Allows datnetW to to be stored with only nrow = nobs(Odata), while datnetA has nrow = p*nobs(Odata) under g_star
      # datnetW <- self$datnetW 
      # datnetA <- self$datnetA
      # sel.sW <- TRUE
      # sel.sA <- TRUE      
      # df.bin.sVar <- self$datnetA$df.bin.sVar
      # print(class(df.bin.sVar)); print(dim(df.bin.sVar)); print(head(df.bin.sVar))
      # print(class(datnetA$df.sVar)); print(dim(datnetA$df.sVar))
      # print(class(datnetW$df.sVar)); print(dim(datnetW$df.sVar))
      # selrow1 <- selrow2 <- selrow3 <- rowsubset
      # if (length(rowsubset) > nrow(df.bin.sVar))     selrow1 <- rowsubset[1:nrow(df.bin.sVar)]
      # if (length(rowsubset) > nrow(datnetA$df.sVar)) selrow2 <- rowsubset[1:nrow(datnetA$df.sVar)]
      # if (length(rowsubset) > nrow(datnetW$df.sVar)) selrow3 <- rowsubset[1:nrow(datnetW$df.sVar)]
      # dfsel <- cbind(df.bin.sVar[selrow1, sel.binsA, drop = FALSE],
      #                 datnetA$df.sVar[selrow2, sel.sA, drop = FALSE],
      #                 datnetW$df.sVar[selrow3, sel.sW, drop = FALSE])
      # VERSION 2: DOESN'T WORK when df.bin.sVar is empty (data.frame()) -> can't cbind 0 rows with N>0 rows
      # Select the columns, cbind together, then select the rows once:
      # dfsel <- cbind(df.bin.sVar[, sel.binsA, drop = FALSE], 
      #                 datnetA$df.sVar[, sel.sA, drop = FALSE], 
      #                 datnetW$df.sVar[, sel.sW, drop = FALSE])
      # print("dfsel: "); print(dim(dfsel))
      # dfsel <- dfsel[rowsubset, , drop = FALSE]
      return(dfsel)
    },
    get.outvar = function(rowsubset = TRUE, var) {
      if (var %in% self$names.sWsA) {
        self$df.sW.sA[rowsubset, var]
      } else if (var %in% colnames(self$df.bin.sVar)) {
        self$df.bin.sVar[rowsubset, var]
      } else {
        stop("outcome variable not found")
      }
    }, 
    # NEW 06/16/15: No longer need to cbind both datasets until last moment. 
    # Can be done with cbind(data.frame(), data.frame()) inside DatNet.sW.sA
    # This function should just return DatNet.sWsA with datnetAstar that cotanis the final matrix dat.sA with n*p rows
    # NEW 06/18/15: It appears that cbind(df1,df2) could be significantly slower and much more error prone.... Keeping both approaches for now.
    # GENERATE SAMPLES OF (dat.sW, dat.sA), where sA is sampled under f.g_name, of size p*nobs, for p>=1.
    # Returns observed (sA,sW) data if is.null(f.g_name)
    # TO ADD: pass ahead a total number of sA that will be created by DatNet class (need it to pre-allocate self$dat.sWsA)
    # TO ADD: Current structure requires building sA twice, once for observed data and once for g_0 when g_0 unknown. This can be expensive. 
    make.df.sWsA = function(p = 1, Kmax, nodes, f.g_name = NULL, f.g_args = NULL)  {
      datnetW <- self$datnetW
      datnetA <- self$datnetA
      assert_that(is.count(p)) # self$p <- p
      nobs <- datnetW$nOdata
      # Copy variable detected types and bin interval definitions from the observed data classes (datnetW, datnetA) to itself:
      # Can assigns types and intervals to itself, since DatNet.sWsA inherits from DatNet
      self$type.sVar <- c(datnetW$type.sVar, datnetA$type.sVar)
      self$cbin_intrvls <- c(datnetW$cbin_intrvls, datnetA$cbin_intrvls)
      if (is.null(f.g_name)) {  # set df.sW.sA to observed data (sW,sA) if g.fun is.null      	
        df.sWsA <- cbind(datnetW$df.sVar, datnetA$df.sVar) # assigning summary measures as data.frames:
      } else {  # sample A from g function in f.g_name
        df.sWsA <- matrix(nrow = (nobs * p), ncol = (datnetW$ncols.sVar + datnetA$ncols.sVar))  # pre-allocate result matx sW.sA
        colnames(df.sWsA) <- self$names.sWsA
        df.sWsA <- data.frame(df.sWsA)
        for (i in seq_len(p)) {
          Avec.df <- data.frame(Anode = f.gen.A.star(Kmax, datnetW$dat.netVar, f.g_name, f.g_args), stringsAsFactors = FALSE)
          colnames(Avec.df)[1] <- nodes$Anode
          # CALL A fun DatNet THAT will instead grow dat.sA in DatNet. # PUT A FLAG IN $new? NOT POSSIBLE. JUST Acess another function on DatNet
          # if (i == 1) {
            datnetAstar <- DatNet$new(Odata = Avec.df, NetInd_k = datnetW$NetInd_k, Kmax = Kmax, nodes = nodes, VarNodes = nodes$Anode)
          # } else {
          #   datnetAstar$add.dat.netVar(Odata = Avec.df)
          # }
          datnetAstar$make_sVar(names.sVar = datnetA$names.sVar) # create summary measures sA
          # Assigin summary measures as data.frames:
          df.sWsA[((i - 1) * nobs + 1):(nobs * i), ] <- cbind(datnetW$df.sVar, datnetAstar$df.sVar)[,]
        }
      }
      self$df.sVar <- df.sWsA
      self$df.sW.sA <- self$df.sVar # for compatibility with older functions
      return(invisible(self))
    }
  ),
  active = list(
    names.sWsA = function() {c(self$datnetW$names.sVar, self$datnetA$names.sVar)},
    nobs = function() { nrow(self$df.sW.sA) }
  ),
  private = list(
    placeholder = list()
  )
)