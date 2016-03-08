#----------------------------------------------------------------------------------
# Class for defining, parsing and evaluating the summary measures.
# Expressions sVar.exprs are evaluated in the environment of a given data.frame.
#----------------------------------------------------------------------------------

is.DefineSummariesClass <- function(obj) "DefineSummariesClass" %in% class(obj)
# Useful function for testing if a name is a valid R object name:
isValidAndUnreservedName <- function(string) {
  make.names(string) == string 	# make.names converts any string into a valid R object name
}

capture.exprs <- function(...) {
  sVar.exprs <- eval(substitute(alist(...)))
  if (length(sVar.exprs)>0) { # deparse into characters when expr is.call, but keep as-is otherwise
    sVar.exprs <- lapply(sVar.exprs, function(x) if (is.character(x)) {x} else {deparse(x)})
  }
  # if not a single argument was named, names attribute of sVar.exprs will be null => add names attribute
  if (is.null(names(sVar.exprs))) names(sVar.exprs) <- rep_len("", length(sVar.exprs))
  if (length(sVar.exprs)!=0 && any(names(sVar.exprs)%in%"")) {
    message("Some summary measures were not named, automatic column name(s) will be generated during evaluation")
  }
  return(sVar.exprs)
}

#' Define Summary Measures sA and sW
#'
#' Define and store summary measures \code{sW} and \code{sA} that can be later processed inside 
#'  \code{\link{eval.summaries}} or \code{\link{tmlenet}} functions.
#'  \code{def_sW} and \code{def_sA} return an \code{R6} object of class \code{\link{DefineSummariesClass}}
#'  which stores the user-defined summary measure functions of the baseline covariates \code{W}
#'  and exposure \code{A}, which can be later evaluated inside the environment of the input \code{data} data frame.
#'  Note that calls to \code{def_sW} must be used for defining the summary measures that are functions 
#'  of \strong{only the baseline covariates \code{W}}, while calls to \code{def_sA} must be used
#'  for defining the summary measures that are functions of both, \strong{the baseline covariates \code{W}
#'  and exposure \code{A}}.
#'  Each summary measure is specified as an evaluable R expression or a string that can be parsed into 
#'  an evaluable R expression. Any variable name that exists as a named column in the input \code{data}
#'  data frame can be used as part of these expressions.
#'  Separate calls to \code{def_sW/def_sA} functions can be aggregated into a single collection with '+' function, 
#'  e.g., \code{def_sW(W1)+def_sW(W2)}.
#'  A special syntax is allowed inside these summary expressions:  
#'  \itemize{
#'  \item \code{'Var[[index]]'} - will index the friend covariate values of the variable \code{Var}, e.g., 
#'    \code{'Var[[1]]'} will pull the covariate value of \code{Var} for the first friend, \code{'Var[[Kmax]]'}
#'    of the last friend, and
#'    \code{'Var[[0]]'} is equivalent to writing \code{'Var'} itself (indexes itself).
#'  }
#'  A special argument named \code{replaceNAw0} can be also passed to the \code{def_sW}, \code{def_sA} functions:
#'  \itemize{
#'  \item \code{replaceNAw0 = TRUE} - automatically replaces all the missing network covariate values 
#'  (\code{NA}) with \code{0}.
#'  }
#'  One can then test the evaluation of these summary measures by either passing the returned 
#'  \code{\link{DefineSummariesClass}} object to function \code{\link{eval.summaries}} or by calling the 
#'  internal method \code{eval.nodeforms(data.df, netind_cl)} on the result returned by \code{def_sW} or \code{def_sA}.
#'  Each separate argument to \code{def_sW} or \code{def_sA} represents a new summary measure.
#'  The user-specified argument name defines the name of the corresponding summary measure 
#'  (where the summary measure represents the result of the evaluation of the corresponding R expression specified by the argument).
#'  When a particular argument is unnamed, the summary measure name
#'  will be generated automatically (see Details, Naming Conventions and Examples below).
#' 
#' @section Details:
#' 
#' The R expressions passed to these functions are evaluated later inside \code{\link{tmlenet}} or 
#'  \code{\link{eval.summaries}} functions,
#'  using the environment of the input data frame, which is enclosed within the user-calling environment.
#' 
#' Note that when observation \code{i} has only \code{j-1} friends, the \code{i}'s value of \code{"W_netFj"} is
#'  automatically set to \code{NA}. 
#'  This can be an undersirable behavior in some circumstances, in which case one can automatically replace all such
#'  \code{NA}'s with \code{0}'s by setting the argument \code{replaceMisVal0 = TRUE} when calling functions 
#'  \code{def_sW} or \code{def_sA}, i.e., \code{def_sW(W[[1]], replaceMisVal0 = TRUE)}.
#' 
#' @section Naming conventions:
#' Naming conventions for summary measures with no user-supplied name (e.g., \code{def_sW(W1)}).
#' 
#' ....................................
#'  \itemize{
#'  \item If only one unique variable name is used in the summary expression (only one parent), use the variable 
#'    name itself to name the summary measure;
#'  \item If there is more than 1 unique variable name (e.g., \code{"W1+W2"}) in the summary expression, throw an exception 
#'    (user must always supply summary measure names for such expressions).
#'  }
#' 
#' Naming conventions for the evaluation results of summary measures defined by \code{def_sW} & \code{def_sA}.
#' 
#' ....................................
#'  \itemize{
#'  \item When summary expression evaluates to a vector result, the vector is first converted to a 1 col matrix,
#'    with column name set equal to the summary expression name;
#'  \item When the summary measure evaluates to a matrix result and the expression has only one unique variable 
#'    name (one parent), the matrix column names are generated as follows: for the expressions such as \code{"Var"}
#'    or \code{"Var[[0]]"}, the column names \code{"Var"} are assigned
#'    and for the expressions such as \code{"Var[[j]]"}, the column names \code{"Var_netFj"} are assigned.
#'  \item When the summary measure (e.g., named \code{"SummName"}) evaluates to a matrix and either: 1) there is
#'    more than one unique variable name used inside the expression (e.g., \code{"A + 2*W"}),
#'    or 2) the resulting matrix has empty (\code{""}) column names, the column names are assigned according to the
#'    convention:
#'    \code{"SummName.1"}, ..., \code{"SummName.ncol"},
#'    where \code{"SummName"} is replaced by the actual summary measure name and \code{ncol} is the number of columns
#'    in the resulting matrix.
#'  }
#' @param ... Named R expressions or character strings that specify the formula for creating the summary measures.
#' @return R6 object of class \code{DefineSummariesClass} which can be passed as an argument to \code{eval.summaries}
#'  and \code{tmlenet} functions.
#' @seealso \code{\link{eval.summaries}} for
#'  evaluation and validation of the summary measures,
#'  \code{\link{tmlenet}} for estimation,
#'  \code{\link{DefineSummariesClass}} for details on how the summary measures are stored and evaluated.
#' @example tests/examples/2_defsWsA_examples.R
#' @export
def_sW <- function(...) {
  # call outside fun that parses ... and assigns empty names "" if the names attribute not set:
  sVar.exprs <- capture.exprs(...)
  sVar.exprs <- c(sVar.exprs, list(nF = "nF")) # add nF node (vector with counts of friends):
  node_evaluator <- DefineSummariesClass$new(type = "sW")
  node_evaluator$set.user.env(user.env = parent.frame())
  node_evaluator$set.new.exprs(exprs_list = sVar.exprs)
  return(node_evaluator)
}

#' @rdname def_sW
#' @export
def_sA <- function(...) {
  # call outside fun that parses ... and assigns empty names "" if the names attribute not set:
  sVar.exprs <- capture.exprs(...)
  node_evaluator <- DefineSummariesClass$new(type = "sA")
  node_evaluator$set.user.env(user.env = parent.frame())
  node_evaluator$set.new.exprs(exprs_list = sVar.exprs)
  return(node_evaluator)
}

#' @rdname def_sW
#' @export
def_new_sA <- function(...) {
  # call outside fun that parses ... and assigns empty names "" if the names attribute not set:
  sVar.exprs <- capture.exprs(...)
  node_evaluator <- DefineSummariesClass$new(type = "new.sA")
  node_evaluator$set.user.env(user.env = parent.frame())
  node_evaluator$set.new.exprs(exprs_list = sVar.exprs)
  return(node_evaluator)
}

# S3 method '+' for adding two DefineSummariesClass objects
# Summary measure lists in both get added as c(,) into the summary measures in sVar1 object
#' @rdname def_sW
#' @param sVar1 An object returned by a call to \code{def_sW} or \code{def_sA} functions.
#' @param sVar2 An object returned by a call to \code{def_sW} or \code{def_sA} functions.
#' @export
`+.DefineSummariesClass` <- function(sVar1, sVar2) {
  assert_that(is.DefineSummariesClass(sVar1))
  assert_that(is.DefineSummariesClass(sVar2))
  assert_that(all.equal(sVar1$type, sVar2$type))
  # remove duplicate nF node from sVar1 (keep the one in sVar2)
  if (sVar1$type %in% "sW") {
    sVar1 <- sVar1$remove.expr(SummaryName = "nF")
  }
  sVar1$add.new.exprs(NewSummaries = sVar2)
  return(sVar1)
}

# ------------------------------------------------------------------------------------------
# Standardize all names (and fill-in the empty names) according TO THE *SAME* *NAMING* *CONVENTION*;
# ------------------------------------------------------------------------------------------
eval.standardize.expr <- function(expr.idx, self, data.df) {
  # browser()
  # -------------------------------------------------------
  # First evaluate the expression result:
  # -------------------------------------------------------
  # eval.nodeform.out_time <- system.time(
    evalres <- eval.nodeform.out(expr.idx = expr.idx, self = self, data.df = data.df)
    # )
  # print("expr_char: " %+% self$exprs_list[[expr.idx]])
  # print("eval.nodeform.out_time: "); print(eval.nodeform.out_time)

  expr_char <- self$exprs_list[[expr.idx]] # expression itself as string
  expr_nm <- names(self$exprs_list)[expr.idx] # current expression name
  expr_parents <- evalres[["par.nodes"]] # names of parents vars for this expression
  # -------------------------------------------------------
  # no user-supplied argument name, hence need to name this expression:
  # -------------------------------------------------------
  # flag TRUE if user did not provide an argument name for self$exprs_list[expr.idx]:
  expr_noname <- (names(self$exprs_list)[expr.idx] %in% "")
  if (expr_noname && (length(expr_parents)>1)) {
    stop("must name complex expressions that involve more than one variable: " %+% expr_char)
  } else if (expr_noname && !is.null(expr_parents) && is.character(expr_parents)) {
    if (gvars$verbose) {
      message("assigning a name '" %+% expr_parents %+% "' to expression: " %+% expr_char)
    }
    expr_nm <- expr_parents
  } else if (expr_noname) stop(expr_char%+% ": parents are null or not a character vector")
  # -------------------------------------------------------
  # convert vectors to columns, name the matrix columns according to the same naming convention:
  # -------------------------------------------------------
  # evaluation result:
  # if result a vector: convert to one-col matrix assign a name: names(self$exprs_list)[expr.idx] = expr_parents
  if (is.vector(evalres[["evaled_expr"]])) {
     expr_res <- matrix(data = evalres[["evaled_expr"]], ncol = 1)
     colnames(expr_res) <- expr_nm
     return(list(new_expr_name = expr_nm, evaled_expr = expr_res, par.nodes = evalres[["par.nodes"]]))
  # for matrix results: if column names exist (!is.null(colnames(expr_res))): DO NOTHING
  # if column names don't exist (is.null(colnames(expr_res))) or some are empty strings "":
  } else if (is.matrix(evalres[["evaled_expr"]])) {
    if (is.null(colnames(evalres[["evaled_expr"]])) || (any(colnames(evalres[["evaled_expr"]])%in%"")) || (length(expr_parents)>1)) {
      # assign names by convention: <- expr_nm%+%"."%+%c(1:ncol(expr_res))
      colnames(evalres[["evaled_expr"]]) <- expr_nm%+%"."%+%c(1:ncol(evalres[["evaled_expr"]]))
    }
    return(list(new_expr_name = expr_nm, evaled_expr = evalres[["evaled_expr"]], par.nodes = evalres[["par.nodes"]]))
  } else {
    # if result is not a vector or matrix: throw an exception
    stop(expr_char%+% ": summary measure result type is " %+%class(evalres[["evaled_expr"]])%+%"; only matrix or a vector results are supported")
  }
}

# -----------------------------------------------------------------------------------------------------------------------
# Special functions for network/time subsetting. 
# Kept in a list, these functions over-ride the standard function `[` and `[[` when evaluating the summary measures
# -----------------------------------------------------------------------------------------------------------------------
node_fun <- list(
  vecapply = function(X, idx, func) { # custom wrapper for apply that turns a vector X into one column matrix
    if (is.vector(X)) dim(X) <- c(length(X), 1) # returns TRUE only if the object is a vector with no attributes apart from names
    # if (is.atomic(x) || is.list(x)) dim(X) <- c(length(X), 1) # alternative way to test for vectors
      x <- parse(text = deparse(func))[[1]]
      nargs <- length(x[[2]])
      if (nargs>1) {
        funline <- deparse(func)[1]
        stop(funline%+%
        ". Node formulas cannot call non-vectorized functions with more than one named argument. If this is a vectorized function, pass its name to set.DAG(, vecfun=).")
      }
    apply(X, idx, func)
  },

  cbind_mod = function(...) { # cbind wrapper for c(,) calls in node formulas, turns one row matrix into repeat Nsamp row matrix
    env <- parent.frame()
    cbind_res <- do.call("cbind", eval(substitute(alist(...)), envir = env) , envir = env)
    if (nrow(cbind_res)==1) {
      # Nsamp <- get("Nsamp", envir = env)
      Nsamp <- env$self$Nsamp
      dprint("env$self$Nsamp:"); dprint(env$self$Nsamp)
      assert_that(!is.null(Nsamp))

      cbind_res <- matrix(cbind_res, nrow = Nsamp, ncol = ncol(cbind_res), byrow = TRUE)
    }
    dprint("cbind_res"); dprint(cbind_res)
    cbind_res
  },

  # custom function for vector look up '['
  # function takes the name of the TD var and index vector => creates a vector of time-varying column names in df
  # returns matrix TD_var[indx]
  # ***NOTE: current '[' cannot evalute subsetting that is based on values of other covariates such as A1C[ifelse(BMI<5, 1, 2)]
  `[` = function(var, indx, ...) {
    env <- parent.frame()
    t <- env$t # t <- get("t", envir = env)
    var <- substitute(var)
    var.chr <- as.character(var)

    if (missing(indx)) stop("missing tindex when using Var[tindex] inside the node formula")
    if (identical(class(indx),"logical")) indx <- which(indx)
    if (is.null(t)) stop("references, s.a. Var[t] are not allowed when t is undefined")
    if (max(indx)>t) stop(paste0(var, "[", max(indx),"] cannot be referenced in node formulas at t = ", t))  # check indx<= t

    # ******* NOTE *******
    # Don't like the current implementation that defines TDvars as characters and then returns a matrix by cbinding 
    # the existing columins in existing data.frame. This is possibly wasteful. Could we instead subset the existing data.frame?
    TDvars <- var.chr%+%"_"%+%indx
    # Checking the variables paste0(var, "_", indx) exist in simulated data.frame environment:
    dprint("ANCHOR_ALLVARNMS_VECTOR_0:"); dprint(env[["ANCHOR_ALLVARNMS_VECTOR_0"]])

    # TO DO: ****
    # EXTEND TO CHECKING FOR TDvar IN ENCLOSING ENVIRONMENT (user.env) AS WELL IF TDvar_t doesn't exist in the data
    # IF TDvar exists check that its a vector of appropriate length, index it accordinly (using which(t%in%tvec))
    # will need to first eval such vector the variable as in:
    # var.val <- eval(var, envir = env)
    existsTDVar <- function(TDvar_t) TDvar_t %in% env[["ANCHOR_ALLVARNMS_VECTOR_0"]]
    check_exist <- sapply(TDvars, existsTDVar)
    if (!all(check_exist)) stop("undefined time-dependent variable(s): "%+%TDvars[which(!check_exist)])
    # THIS STEP COULD BE MORE MEMORY EFFICIENT IF WAS SUBSETTING INSTEAD (BY COLS) ON EXISTING data MATRIX:
    TDvars_eval <- eval(parse(text=paste0("cbind(",paste0(TDvars, collapse=","),")")), envir = env)
    return(TDvars_eval)
  },

  # Builds netVar matrix by using matrix env$NetIndobj$NetInd_k, cbind on result
  # For W[[0]] to work without if else below need to do this:
  # NetInd_k <- cbind(c(1:n), NetInd_k) and then netidx <- netidx + 1
  `[[` = function(var, netidx, ...) {
    env <- parent.frame()
    t <- env$t # t <- get("t", envir = env)
    if (!is.null(t)) stop("simultaneous time varying node references Var[t] and network references Var[[netidx]] are currently not supported")
    if (missing(netidx)) stop("network index (netidx) must be specified when using Var[[netidx]]")
    netind_cl <- env$netind_cl
    if (is.null(netind_cl)) stop("Network must be defined when using Var[[netidx]] syntax")
    Kmax <- netind_cl$Kmax

    var <- substitute(var)
    var.chr <- as.character(var)
    if (! (var.chr %in% env[["ANCHOR_ALLVARNMS_VECTOR_0"]])) stop("variable " %+% var.chr %+% " doesn't exist")
    var.val <- eval(var, envir = env)
    n <- length(var.val)
    if (identical(class(netidx),"logical")) netidx <- which(netidx)
    netVars_eval <- matrix(0L, nrow = n, ncol = length(netidx))
    colnames(netVars_eval) <- netvar(var.chr, netidx)

    # make_net_mat_time <- system.time(
    for (neti in seq_along(netidx)) {
      if (netidx[neti] %in% 0L) {
        netVars_eval[, neti] <- var.val
      } else {
        netVars_eval[, neti] <- var.val[netind_cl$NetInd_k[, netidx[neti]]]
        # opting for replace on entire netVars_eval, will need to do benchmarks later to compare:
        # netVars_eval[is.na(netVars_eval[, neti]), neti] <- env$misXreplace
      }
    }
    # )
    # print("make_net_mat_time: "); print(make_net_mat_time)

    # Don't need to do this if env$misXreplace==gvars$misval (i.e., when want to leave NAs as is)
    # replaceNA_time <- system.time(
      netVars_eval[is.na(netVars_eval)] <- env$misXreplace
      # )
    # print("replaceNA_time"); print(replaceNA_time)
    
    return(netVars_eval)
  }
)

update.intervention.sA <- function(new.sA, sA) {
  assert_that(is.DefineSummariesClass(new.sA))
  # ------------------------------------------------------------------------------------------------
  # DEFINE AND EVALUATE INTERVENTION SUMMARIES:
  # ------------------------------------------------------------------------------------------------
  # 1. Copy (clone) sA object since we are going to modify it (over-write some summaries in sA$exprs_list with new intervention summaries)
  intervened.sA <- sA$clone()
  intervened.sA$type <- new.sA$type
 
  # 2. Replace summaries/nodes in observed sA that were also defined in new.sA with their new (intervened) expressions:
  for (expr.name.sA in names(new.sA$exprs_list)) {
    intervened.sA$replace.expr(SummaryName = expr.name.sA, newSummaries = new.sA)
  }
  new.sA <- intervened.sA
  print("updated expression list for new.sA: "); print(new.sA$exprs_list)
  Anodes <- new.sA$Anodes

  print("detected Anodes based on new.sA:"); print(Anodes)

  # 3. Evaluate a copy of new summaries ->
    # Note that this will require making a copy of OdataDT in OdataDT_R6, since calling
    # datnetA$make.sVar(Odata = OdataDT_R6, sVar.object = new.sA) will modify OdataDT IN PLACE
    # in all instances at the same time: datnetW$dat.sVar, datnetA$dat.sVar & DatNet.ObsP0$dat.sVar
    # OdataDT_R6$OdataDT
    # OdataDT_R6$backupAnodes(Anodes = Anodes)
    # OdataDT_R6$A_g0_DT
    # new.sA.dat <- new.sA$eval.nodeforms(data.df = OdataDT_R6$OdataDT, netind_cl = netind_cl)
    
    # OdataDT_R6$restoreAnodes(Anodes = Anodes)
    # OdataDT_R6$OdataDT
    # obs.sA.dat <- sA$eval.nodeforms(data.df = OdataDT_R6$OdataDT, netind_cl = netind_cl)

  # 4. Clean-up input data: remove all variables not mentioned in sW/sA from OdataDT_R6$OdataDT
    # This could be a problem if a formula, such as def_sw(PA0 = (PA==0)) exists, but PA itself wasn't defined
    # The code below will then remove PA from the observed data (which should have been kept)
    # *** NEED TO ALSO ADD Y TO THIS (outcome)
    # allvars <- as.vector(c(unlist(sW$sVar.names.map), unlist(sA$sVar.names.map), "Y"))
    # allDatNames <- names(OdataDT_R6$OdataDT)
    # remove_vars <- allDatNames[!(allDatNames %in% allvars)]
    # OdataDT_R6$OdataDT[, (remove_vars):=NULL, with = FALSE]
  return(new.sA)
}

## ---------------------------------------------------------------------
#' R6 class for parsing and evaluating user-specified summary measures (in \code{exprs_list})
#'
#' This \pkg{R6} class that inherits from \code{Define_sVar} and can parse and evaluate (given the input data frame) the summary measures defined by functions 
#'  \code{\link{def_sW}} and \code{\link{def_sA}}. 
#'  The object of this class is generally instantiated by calling functions \code{def_sA} or \code{def_sW}.
#'  The summary expressions (stored in \code{exprs_list}) are evaluated in the environment of the input data.frame.
#'  Note that the evaluation results of the summary measures are never stored inside this class, 
#'  data can be stored only inside \code{\link{DatNet}} and \code{\link{DatNet.sWsA}} \pkg{R6} classes.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{type}} - Type of the summary measure, \code{sW} or \code{sA}, determined by the calling functions \code{\link{def_sW}} or \code{\link{def_sA}}.
#' \item{\code{exprs_list}} - Deparsed list of summary expressions (as strings).
#' \item{\code{new_expr_names}} - The summary measure names, if none were provided by the user these will be 
#'  evaluated on the basis of variable names used in the summary expression itself.
#' \item{\code{sVar.names.map}} - Named list that maps the user specified summary measure names to the corresponding matrix column names
#'  from the summary measure evaluation result.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(type)}}{Instantiate a new object of class \code{DefineSummariesClass} by providing a type, \code{"sW"} or \code{"sA"}.}
#'   \item{\code{set.new.exprs(exprs_list)}}{Sets the internal summary measure expressions to the list provided in \code{exprs_list}.}
#'   \item{\code{add.new.exprs(NewSummaries)}}{Adds new internal summary measure expressions to the existing ones, \code{NewSummaries} 
#'    must be an object of class \code{DefineSummariesClass} (to enable \code{Object1 + Object2} syntax).}
#'   item{\code{remove.expr(SummaryName)}}{Remove expression by name (for removing duplicate 'nF' expressions for repeated calls with def_sW()+def_sW()).}
#'   \item{\code{eval.nodeforms(data.df, netind_cl)}}{Evaluate the expressions one by one, standardize all names according to one naming
#'    convention (described in \code{\link{def_sW}}), \code{cbind}ing results together into one output matrix. \code{data.df} is the input 
#'    data.frame and \code{netind_cl} is the input network stored in an object of class \code{\link[simcausal]{NetIndClass}}.}
#'   \item{\code{df.names(data.df)}}{List of variables in the input data \code{data.df} gets assigned to a special
#'    variable (\code{ANCHOR_ALLVARNMS_VECTOR_0}).}
#' }
#' @importFrom assertthat assert_that
#' @export
DefineSummariesClass <- R6Class("DefineSummariesClass",
  class = TRUE,
  portable = TRUE,
  # inherit = Define_sVar,
  public = list(
    type = NA,                    # "sW" or "sA" depending on which functions called the constructor: def_sW or def_sA
    exprs_list = list(),          # list of expressions in character strings, with attribute "names" set to the user-supplied names of the expressions

    user.env = NULL,              # user environment to be used as enclos arg to eval(sVar)
    cur.node = list(),            # current evaluation node (set by self$eval.nodeforms())
    asis.flags = list(),          # list of flags, TRUE for "as is" node expression evaluation
    ReplMisVal0 = FALSE,          # vector of indicators, for each TRUE sVar.expr[[idx]] will replace all NAs with gvars$misXreplace (0)
    sVar.misXreplace = NULL,      # replacement values for missing sVar, vector of length(exprs_list)
    netind_cl = NULL,
    Kmax = NULL,
    Nsamp = NULL,                 # sample size (nrows) of the simulation dataset

    new_expr_names = list(),      # re-evaluated summary measure names, if non provided by the user these will be evaluated on the basis of variable names used in the expression
    sVar.names.map = list(),      # the map between user-supplied expression (argument names) to the column names of each expression in self$exprs_list

    Anodes = character(),         # vector with unique intervention node/summary names (from new.sA)

    node_fun = node_fun,
    evaled_before = FALSE,        # Indicator that these summaries have been evaluated on the real data at least once

    initialize = function(type) {
      self$type <- type
      invisible(self)
    },

    # define new summary measures to be evaluated later:
    # will define 1) self$exprs_list; 2) self$asis.flags; 3) self$ReplMisVal0; 4) self$sVar.misXreplace
    # initialize the map from self$exprs_list to variable names in each expression (column names)
    set.new.exprs = function(exprs_list) {
      self$exprs_list <- exprs_list
      self$asis.flags <- attributes(exprs_list)[["asis.flags"]]

      # check for special argument replaceNAw0, if exists, remove it from the list of expressions:
      if (any(names(self$exprs_list) %in% "replaceNAw0")) {
        ReplMisVal0.idx <- which(names(self$exprs_list) %in% "replaceNAw0")
        self$ReplMisVal0 <- as.logical(self$exprs_list[[ReplMisVal0.idx]])
        self$exprs_list <- self$exprs_list[-ReplMisVal0.idx]
        if (gvars$verbose) {
          print("Detected replaceNAw0 flag with value: " %+% self$ReplMisVal0);
        }
      }
      # If doesn't already exist, init setting for the names attribute of self$exprs_list:
      if (is.null(names(self$exprs_list))) names(self$exprs_list) <- rep_len("", length(self$exprs_list))

      # Error when un-named expressions (name(self$exprs_list[[i]])=="") are defined for intervention summaries:
      if ((self$type %in% "new.sA")  && length(self$exprs_list) != 0 && (is.null(names(self$exprs_list)) || any(names(self$exprs_list)==""))) {
        stop("must provide a name for each node expression")
      }

      if (is.null(self$asis.flags)) {
        self$asis.flags <- as.list(rep.int(FALSE, length(self$exprs_list)))
        names(self$asis.flags) <- names(self$exprs_list)
      }

      self$ReplMisVal0 <- rep_len(self$ReplMisVal0, length(self$exprs_list))
      self$sVar.misXreplace <- ifelse(self$ReplMisVal0, gvars$misXreplace, gvars$misval)

      # self$sVar.noname <- rep_len(self$sVar.noname, length(self$exprs_list))
      if (gvars$verbose) {
        print("Final node expression(s) list: "); print(self$exprs_list)
      }
      self$sVar.names.map <- vector(mode="list", length = length(self$exprs_list))
      invisible(self)
    },

    # add summary measures to existing ones (to enable Object1 + Object2 syntax):
    add.new.exprs = function(NewSummaries) {
      assert_that(is.DefineSummariesClass(NewSummaries))
      self$exprs_list <- c(self$exprs_list, NewSummaries$exprs_list)
      self$asis.flags <- c(self$asis.flags, NewSummaries$asis.flags)
      self$ReplMisVal0 <- c(self$ReplMisVal0, NewSummaries$ReplMisVal0)
      self$sVar.misXreplace <- c(self$sVar.misXreplace, NewSummaries$sVar.misXreplace)
      self$sVar.names.map <- c(self$sVar.names.map, NewSummaries$sVar.names.map)
      return(self)
    },

    # remove existing summary measure
    # to enable overwriting nF summary measure in Object1 with nF summary measure from Object2
    # when doing Object1 + Object2 syntax.
    remove.expr = function(SummaryName) {
      assert_that(is.character(SummaryName) && (length(SummaryName)==1L) && (!SummaryName%in%""))
      if (any(names(self$exprs_list) %in% SummaryName)) {
        remove_idx <- which(names(self$exprs_list)%in% SummaryName)
        self$exprs_list <- self$exprs_list[-remove_idx]
        self$asis.flags <- self$asis.flags[-remove_idx]
        self$ReplMisVal0 <- self$ReplMisVal0[-remove_idx]
        self$sVar.misXreplace <- self$sVar.misXreplace[-remove_idx]
        self$sVar.names.map <- self$sVar.names.map[-remove_idx]
      }
      return(self)
    },

    # Replace old summary measure with its new definition (give an error if the summary doesn't exist)
    # to enable overwriting observed summaries/nodes in def_sA with new intervention summaries/nodes.
    replace.expr = function(SummaryName, newSummaries) {
      assert_that(is.character(SummaryName) && (length(SummaryName)==1L) && (!SummaryName%in%""))

      if (!self$evaled_before) stop("Can't replace existing summaries which haven't yet been evaluated on the real data at least once")

      newSummaries_idx <- which(names(newSummaries$exprs_list) %in% SummaryName)

      if (any(names(self$exprs_list) %in% SummaryName)) {

        replace_idx <- which(names(self$exprs_list) %in% SummaryName)
        self$exprs_list[replace_idx] <- newSummaries$exprs_list[newSummaries_idx]
        self$asis.flags[replace_idx] <- newSummaries$asis.flags[newSummaries_idx]
        self$ReplMisVal0[replace_idx] <- newSummaries$ReplMisVal0[newSummaries_idx]
        self$sVar.misXreplace[replace_idx] <- newSummaries$sVar.misXreplace[newSummaries_idx]
        # self$sVar.names.map[replace_idx] <- newSummaries$sVar.names.map[newSummaries_idx] # should not be replacing the map

        # mark the summary/node name that was replaced:
        self$Anodes <- c(self$Anodes, SummaryName)

      } else {
        stop("the intervention summary measure name has not been previously defined: " %+% SummaryName)
      }

      return(self)
    },

    # Evaluate the expressions one by one, standardize all names according to one naming convention,
    # cbinding results together into one output matrix
    eval.nodeforms = function(data.df, netind_cl) {
      assert_that(is.data.frame(data.df) | is.data.table(data.df))
      if (is.data.frame(data.df)) setDT(data.df)

      if (missing(netind_cl) && is.null(self$netind_cl)) stop("must specify netind_cl arg at least once")
      if (!missing(netind_cl)) self$netind_cl <- netind_cl
      self$Kmax <- self$netind_cl$Kmax
      self$Nsamp <- nrow(data.df)

      sVar.res_l <- self$new_expr_names <- self$sVar.names.map <- vector(mode = "list", length = length(self$exprs_list))

      # eval_and_addDT <- system.time(
        for (i in seq_along(self$exprs_list)) {
          sVar.eval.res <- eval.standardize.expr(i, self = self, data.df = data.df)
          # sVar.res_l[[i]] <- sVar.eval.res

          self$new_expr_names[[i]] <- sVar.eval.res$new_expr_name
          self$sVar.names.map[[i]] <- colnames(sVar.eval.res$evaled_expr)

          for (colname in colnames(sVar.eval.res$evaled_expr))
            data.df[, (colname):= sVar.eval.res$evaled_expr[,colname]]
        }
      # )

      # print("eval_and_addDT"); print(eval_and_addDT)
      # eval_only <- system.time(sVar.res_l <- lapply(seq_along(self$exprs_list), eval.standardize.expr, self = self, data.df = data.df))
      # print("eval_only: "); print(eval_only)

      # self$new_expr_names <- lapply(sVar.res_l, '[[', 'new_expr_name')
      names(self$new_expr_names) <- unlist(self$new_expr_names)
      names(self$exprs_list) <- names(self$new_expr_names)
      # names(sVar.res_l) <- names(self$new_expr_names)

      # assign self$sVar.names.map based on newly standardized summary names:
      # self$sVar.names.map <- lapply(sVar.res_l, function(x) colnames(x[["evaled_expr"]]))
      names(self$sVar.names.map) <- names(self$new_expr_names)

      # ************************************************************************************
      # CHANGE THIS TO MERGING DUPLICATE SUMMARY MEASRURES INTO ONE
      # ************************************************************************************
      # 1) remove all duplicate summary measures (by name), keeping the ones that were added last:
      if (length(unique(names(self$new_expr_names))) < length(names(self$new_expr_names))) {
        duplic_idx <- duplicated(self$new_expr_names, fromLast = TRUE)
        message("warning: detected duplicate summary measure names, (" %+%
                paste0(self$new_expr_names[duplic_idx], collapse=",") %+%
                "), all duplicates starting from first to last will be removed...")
        # sVar.res_l <- sVar.res_l[-duplic_idx]
        self$sVar.names.map <- self$sVar.names.map[-duplic_idx]
        self$new_expr_names <- self$new_expr_names[-duplic_idx]
        self$exprs_list <- self$exprs_list[-duplic_idx]
      }
      # mat.sVar <- do.call("cbind", lapply(sVar.res_l, function(x) x[["evaled_expr"]]))

      # Indicator that these summaries have been evaluated on the real data at least once:
      self$evaled_before <- TRUE

      # ************************************************************************************
      # THIS IS NOW DONE AUTOMATICALLY WHEN ASSIGNING WITHIN data.table data.df
      # ************************************************************************************
      # 2) remove duplicate columns, keeping the ones that were added last:
      # if (length(unique(colnames(mat.sVar))) < length(colnames(mat.sVar))) {
      #   duplic_idx <- duplicated(colnames(mat.sVar), fromLast = TRUE)
      #   message("warning: detected duplicate column names in summary evaluation matrix, (" %+%
      #           paste0(colnames(mat.sVar)[duplic_idx], collapse=",") %+%
      #           "), all duplicates starting from first to last will be removed...")
      #   mat.sVar <- mat.sVar[,!duplic_idx]
      # }
      # return(mat.sVar)

      return(invisible(data.df))
    },

    # This user.env is used for eval'ing each sVar exprs (enclos = user.env)
    set.user.env = function(user.env) {
      assert_that(!is.null(user.env))
      assert_that(is.environment(user.env))
      self$user.env <- user.env
    },

    # List of variable names from data.df with special var name (ANCHOR_ALLVARNMS_VECTOR_0):
    df.names = function(data.df) {
      return(list(ANCHOR_ALLVARNMS_VECTOR_0 = colnames(data.df)))
    }
  ),

  active = list(
    placeholder = function() {}
  ),

  private = list(
    privplaceholder = function() {}
  )
)
