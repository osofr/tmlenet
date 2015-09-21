
#----------------------------------------------------------------------------------
# Class for defining, parsing and evaluating the summary measures.
# Expressions sVar.exprs are evaluated in the environment of a given data.frame.
#----------------------------------------------------------------------------------
  # II) sVar naming:
  #todo 34 (sVar_evaluator, sVar.name) +0: Instead of throwing an error for un-named vector results (e.g., def.sW(A)), 
    # would be nice to extract the first is.name of the expression and use that as a name.
  #todo 7 (sVar_evaluator, sVar.name) +0: Check that the resulting column names in sVar are all unique!
  # Create an option $keep.sVar.nms; When TRUE do not change the output column names for sVar mat with ncol > 1. 
  # Create an option to overwrite sVar mat colname(s) with user-provided names
  #todo 28 (sVar_evaluator) +0: Consider returning sVar.res_l instead of mat.sVar, also see if there are faster alternatives to cbind
    # (i.e., pre-allocating sVar.mat); perform benchmarks to see if there is any noticable benefit
  #todo 42 ('+.Define_sVar') +0: Allow adding character vector summary measures for sVar2, s.a., 
  # def.sW(W2[[1:Kmax]]) + "netW3_sum = rowSums(W3[[1:Kmax]]"

# Useful function for testing if a name is a valid R object name:
isValidAndUnreservedName <- function(string) {
  make.names(string) == string 	# make.names converts any string into a valid R object name
}

# Wrappers for Define_sVar$new(...) constructor:
#' Define Summary Measures sA and sW
#'
#' Define and store \code{tmlenet} function arguments \code{sW} and \code{sA},
#'  which are arbitrary summary measures based on treatmet \code{Anode} and baseline covariates in \code{Wnodes} 
#'  constructed from the input \code{data}.
#'  Each summary measure \code{sVar} is specified as a named R expression or a named character argument. 
#'  Separate calls to \code{def.sW/def.sA} functions can be aggregated into a single collection with '+' function, 
#'  e.g., \code{def.sW(W1=W1)+def.sW(W2=W2)}.
#'  Additional named arguments that can be passed to \code{def.sW/def.sA} are:
#'  \itemize{
#'  \item \code{noname = TRUE} - to not use the summary measure name provided by the user when assigning
#'    the names to the summary measure columns (variables); and
#'  \item \code{replaceNAw0 = TRUE} - To automatically replace all the missing network covariate values (\code{NA}s) with \code{0}.
#' }
#' @section Details: 
#' The R expressions passed to these functions are evaluated later inside \code{\link{tmlenet}} or \code{\link{eval.summaries}} functions,
#' using the environment of the input data frame and enclosed within the user-calling environment.
#' @param ... Named R expressions or character strings that specify the formula for creating the summary measures.
#' @return R6 object of class \code{Define_sVar} which must be passed as argument to \code{\link{tmlenet}}.
#' @seealso \code{\link{tmlenet}}
#' @example tests/examples/2_defsWsA_examples.R
#' @export
def.sW <- function(...) Define_sVar$new(..., type = "sW", user.env = parent.frame())
#' @rdname def.sW
#' @export
def.sA <- function(...) Define_sVar$new(..., type = "sA", user.env = parent.frame())

is.Define_sVar <- function(obj) "Define_sVar" %in% class(obj)

# S3 method '+' for adding two Define_sVar objects
# Summary measure lists in both get added as c(,) into the summary measures in sVar1 object
#' @rdname def.sW
#' @param sVar1 An object returned by a call to \code{def.sW} or \code{def.sA} functions.
#' @param sVar2 An object returned by a call to \code{def.sW} or \code{def.sA} functions.
#' @export
`+.Define_sVar` <- function(sVar1, sVar2) {
  assert_that(is.Define_sVar(sVar1))
  assert_that(is.Define_sVar(sVar2))
  assert_that(sVar1$type == sVar2$type)
  # if (!is.Define_sVar(sVar1) && is.Define_sVar(sVar2)) {
  #   tmp <- obj1
  #   obj1 <- obj2
  #   obj2 <- tmp
  # }
  # if ("DAG.action" %in% class(obj2)) {
    ## Adding action, possibly with attributes
    ## Option 1: Non-named argument defines the action name
    ## Option 2: Non-named argument defines the nodes
    # if (is.null(names(obj2))) {
    #   toadd <- unlist(obj2, recursive=FALSE)
    #   attr <- list()
    # } else {
      # toadd <- unlist(obj2[names(obj2)==""])
    # name <- unlist(obj2[names(obj2)==""], recursive=FALSE)
    # if (length(name)>1) stop("only one unnamed argument can be specified")
    # if (length(name)==0) name <- obj2[[which(names(obj2)%in%"name")]]
    # if (length(name)==0) stop("name argument for action must be specified")
    # nodes <- unlist(obj2[names(obj2)%in%"nodes"])
    # nodes <- unlist(obj2[names(obj2)%in%"nodes"], recursive=FALSE)
    # nodes <- obj2[names(obj2)%in%"nodes"][[1]]
    # attr <- obj2[(names(obj2)!="") & (!(names(obj2) %in% c("name", "nodes")))]
    # dprint("name"); dprint(name)
    # dprint("nodes"); dprint(nodes)
    # dprint("attr"); dprint(attr)
    # res <- add.action(DAG=obj1, name=name, nodes=nodes, attr=attr)
  # } else if ("DAG.nodelist" %in% class(obj2)) {
    # res <- c(obj1, obj2)
    # class(res) <- "DAG"
    # res <- add.nodes(DAG=obj1, nodes=obj2)
  # } else {
    # stop("Cannot add unknown type to DAG")
  # }
  sVar1$ReplMisVal0 <- c(sVar1$ReplMisVal0, sVar2$ReplMisVal0)
  sVar1$sVar.misXreplace <- c(sVar1$sVar.misXreplace, sVar2$sVar.misXreplace)
  sVar1$sVar.exprs <- c(sVar1$sVar.exprs, sVar2$sVar.exprs)
  sVar1$sVar.expr.names <- c(sVar1$sVar.expr.names, sVar2$sVar.expr.names)
  sVar1$sVar.noname <- c(sVar1$sVar.noname, sVar2$sVar.noname)
  # res <- add.action(DAG=obj1, name=name, nodes=nodes, attr=attr)
  # res
  return(sVar1)
}

# take sVar expression index and evaluate:
parse.sVar.out <- function(sVar.idx, self) {
  sVar.expr <- self$sVar.exprs[[sVar.idx]]
  sVar.name <- self$sVar.expr.names[sVar.idx]
  misXreplace <- self$sVar.misXreplace[sVar.idx]
  # ******
  eval.sVar.params <- c(self$df.names(self$data.df), list(misXreplace = misXreplace), list(netind_cl = self$netind_cl))
  data.env <- c(eval.sVar.params, self$node_fun, self$data.df)
  # ******

  if (is.character(sVar.expr)) {
    sVar.expr_call <- try(parse(text=sVar.expr)[[1]]) 	# parse expression into a call
    if(inherits(sVar.expr_call, "try-error")) {
      stop("error while evaluating expression: " %+% sVar.expr %+% ".\nCheck syntax specification.", call.=FALSE)
    }
  } else if (is.call(sVar.expr)){
    sVar.expr_call <- sVar.expr
    message(sVar.expr_call %+% ": sVar formula is already a parsed call")
  } else {
    stop("sVar formula class: " %+% class(sVar.expr) %+% ". Currently can't process sVar formulas that are not strings or calls.")
  }
  sVar.expr_call <- eval(substitute(substitute(e, list(Kmax = eval(self$Kmax))), list(e = sVar.expr_call))) # Replace Kmax its val
  evaled_expr <- try(eval(sVar.expr_call, envir = data.env, enclos = self$user.env)) # eval'ing expr in the envir of data.df

  no.sVar.name <- self$sVar.noname[sVar.idx]
  # no.sVar.name <- is.null(sVar.name) || (sVar.name %in% "")

  if (is.matrix(evaled_expr)) {
    if (no.sVar.name) sVar.name <- colnames(evaled_expr)
    if (!no.sVar.name && ncol(evaled_expr) > 1) sVar.name <- sVar.name %+% "." %+% (1 : ncol(evaled_expr))
    # if (no.sVar.name || ncol(evaled_expr) > 1) message(sVar.expr %+% ": 
      # the result matrix is assigned the following column name(s): " %+% paste(sVar.name, collapse = ","))
  } else {
    if (no.sVar.name) stop(sVar.expr %+% ": summary measures not defined with Var[[...]] must be named.")
    evaled_expr <- as.matrix(evaled_expr)
  }
  colnames(evaled_expr) <- sVar.name
  return(evaled_expr)
}


## ---------------------------------------------------------------------
#' R6 class for parsing and evaluating user-specified summary measures (in \code{sVar.exprs})
#'
#' This R6 class can parse and evaluate (given the input data) the summary measures defined with functions 
#'  \code{\link{def.sW}} and \code{\link{def.sA}}. 
#'  The summary measure expressions (stored in \code{sVar.exprs}) are evaluated in the environment of the input data.frame.
#'  Note that the resulting summary measures are never stored inside this class, 
#'  the data is stored only in \code{\link{DatNet}} and \code{\link{DatNet.sWsA}} R6 classes.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{data.df}} - Input data frame that will be used for evaluating the summary measure expressions.
#' \item{\code{Kmax}} - Maximum number of friends for any observation.
#' \item{\code{netind_cl}} - Pointer to a network instance of class \code{simcausal::NetIndClass}.
#' \item{\code{ReplMisVal0}} - A logical vector that captures args \code{replaceNAw0=TRUE/FALSE} in functions \code{def.sW}, \code{def.sA}.
#'  If \code{TRUE} for a particular summary measure in \code{sVar.exprs} then all missing network \code{VarNode} 
#'  values (when \code{nF[i] < Kmax}) will get replaced with \code{tmlenet:::gvars$misXreplace} (default is 0).
#' \item{\code{sVar.misXreplace}} - Maximum number of friends for any observation.
#' \item{\code{sVar.noname}} - Vector, for each \code{TRUE}, \code{sVar.expr[[idx]]} ignores user-supplied name and generates names automatically.
#' \item{\code{sVar.exprs}} - Deparsed summary measure expressions expressions (character vector).
#' \item{\code{sVar.expr.names}} - User-provided names for each summary measure in \code{sVar.expr}.
#' \item{\code{sVar.names.map}} - Map of the summary measure names given by the user to its actual column names
#'  (one name can represent a multivariate summary measure).
#' \item{\code{type}} - Type of the summary measure, sW or sA, determined by the calling functions \code{def.sW} or \code{def.sA}.
#' \item{\code{user.env}} - Captured user-environment from which \code{def.sW} or \code{def.sA} was called.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(..., type, user.env)}}{...}
#'   \item{\code{get.mat.sVar(data.df, netind_cl, addnFnode = NULL)}}{...}
#'   \item{\code{df.names(data.df)}}{List of variables in the input data \code{data.df} gets assigned to a special 
#'   variable (\code{ANCHOR_ALLVARNMS_VECTOR_0})}
#' }
#' @importFrom assertthat assert_that
#' @export
Define_sVar <- R6Class("Define_sVar",
  class = TRUE,
  portable = TRUE,
  public = list(
    data.df = NULL,               # data.frame that is used for evaluation of sVar expressions (passed to get.mat.sVar)
    Kmax = NULL,
    netind_cl = NULL,
    ReplMisVal0 = FALSE,          # Replace missing network VarNode values (when nF[i] < Kmax) with gvars$misXreplace (0)?
    sVar.misXreplace = NULL,      # Replacement values for missing sVar (length(sVar.exprs)), either gvars$misXreplace or gvars$misval
    sVar.noname = FALSE,          # Vector, for each TRUE sVar.expr[[idx]] ignores user-supplied name and generates names automatically
    sVar.exprs = character(),     # deparsed sVar expressions (char vector)
    sVar.expr.names = character(),# user-provided name of each sVar.expr
    sVar.names.map = list(),
    type = NA,                    # sW or sA
    user.env = emptyenv(),        # User environment used as enclos arg to eval(sVar, enclos=)
    # mat.sVar = matrix(),        # no longer storing the sVar evaluation result

    node_fun = list(
      # Builds netVar matrix by using matrix env$NetIndobj$NetInd_k, cbind on result
      # For W[[0]] to work without if else below need to do this:
      # NetInd_k <- cbind(c(1:n), NetInd_k) and then netidx <- netidx + 1
      `[[` = function(var, netidx, ...) {
        env <- parent.frame()
        netind_cl <- env$netind_cl
        if (missing(netidx)) stop("network index (netidx) must be specified when using Var[[netidx]]")
        var <- substitute(var)
        var.chr <- as.character(var)
        if (! (var.chr %in% env[["ANCHOR_ALLVARNMS_VECTOR_0"]])) stop("variable " %+% var %+% " doesn't exist")
        var.val <- eval(var, envir = env)
        n <- length(var.val)
        if (identical(class(netidx),"logical")) netidx <- which(netidx)
        netVars_eval <- matrix(0L, nrow = n, ncol = length(netidx))
        colnames(netVars_eval) <- netvar(var.chr, netidx)
        for (neti in seq_along(netidx)) {
          if (netidx[neti] %in% 0L) {
            netVars_eval[, neti] <- var.val
          } else {
            netVars_eval[, neti] <- var.val[netind_cl$NetInd_k[, netidx[neti]]]
            # opting for replace on entire netVars_eval, will need to do benchmarks later to compare:
            # netVars_eval[is.na(netVars_eval[, neti]), neti] <- env$misXreplace
          }
        }
        # Don't need to do this if env$misXreplace==gvars$misval (i.e., when want to leave NAs as is)
        netVars_eval[is.na(netVars_eval)] <- env$misXreplace
        return(netVars_eval)
      }
    ),

    # capture sVar expressions and capture the user environment;
    # user.env is used when eval'ing sVar exprs (enclos = user.env)
    initialize = function(..., type, user.env) {
      self$type <- type
      self$user.env <- user.env
      self$sVar.exprs <- eval(substitute(alist(...)))
      if (length(self$sVar.exprs)>0) { # deparse into characters when expr is.call, but keep as-is otherwise
        self$sVar.exprs <- lapply(self$sVar.exprs, function(x) if (is.character(x)) {x} else {deparse(x)})
      }
      self$sVar.expr.names <- names(self$sVar.exprs)
      if (is.null(self$sVar.expr.names)) self$sVar.expr.names <- rep_len("", length(self$sVar.exprs))
      if (length(self$sVar.exprs) != 0 && (is.null(self$sVar.expr.names) || any(self$sVar.expr.names==""))) {
        stop("must provide a name for each summary measure expression")
        # message("Some summary measures were not named, automatic column name(s) will be generated during evaluation")
      }

      if (any(self$sVar.expr.names %in% "replaceNAw0")) {
        ReplMisVal0.idx <- which(self$sVar.expr.names %in% "replaceNAw0")
        self$ReplMisVal0 <- as.logical(self$sVar.exprs[[ReplMisVal0.idx]])
        self$sVar.expr.names <- self$sVar.expr.names[-ReplMisVal0.idx]
        self$sVar.exprs <- self$sVar.exprs[-ReplMisVal0.idx]
        message("Detected replaceNAw0 flag with value: " %+% self$ReplMisVal0);
      }
      if (any(self$sVar.expr.names %in% "noname")) {
        noname.idx <- which(self$sVar.expr.names %in% "noname")
        self$sVar.noname <- as.logical(self$sVar.exprs[[noname.idx]])
        self$sVar.expr.names <- self$sVar.expr.names[-noname.idx]
        self$sVar.exprs <- self$sVar.exprs[-noname.idx]
        message("Detected noname flag with value: " %+% self$sVar.noname);
      }

      self$ReplMisVal0 <- rep_len(self$ReplMisVal0, length(self$sVar.exprs))
      self$sVar.misXreplace <- ifelse(self$ReplMisVal0, gvars$misXreplace, gvars$misval)
      self$sVar.noname <- rep_len(self$sVar.noname, length(self$sVar.exprs))

      # for later S3 dispatch?
      if (!is.na(type)) class(self) <- c(class(self), type)
      message("Defined the following summary measure(s): "); print(self$sVar.exprs)
      invisible(self)
    },

    get.mat.sVar = function(data.df, netind_cl, addnFnode = NULL) {
      assert_that(is.data.frame(data.df))
      self$data.df <- data.df
      if (missing(netind_cl) && is.null(self$netind_cl)) stop("must specify netind_cl arg at least once")
      if (!missing(netind_cl)) self$netind_cl <- netind_cl
      self$Kmax <- self$netind_cl$Kmax
      # call lapply on parse.sVar.out for each sVar in sVar.expr.names -> sVar.res_l
      sVar.res_l <- lapply(seq_along(self$sVar.exprs), parse.sVar.out, self = self)
      names(sVar.res_l) <- self$sVar.expr.names

      if (!is.null(addnFnode)) sVar.res_l <- c(sVar.res_l, list(nF = netind_cl$mat.nF(addnFnode)))

      # SAVE THE MAP BETWEEEN EXPRESSION NAMES AND CORRESPONDING COLUMN NAMES:
      self$sVar.names.map <- lapply(sVar.res_l, colnames)
      mat.sVar <- do.call("cbind", sVar.res_l)
      self$data.df <- NULL # erase pointer to data after evaluation is done
      return(mat.sVar)
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