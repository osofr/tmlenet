
#----------------------------------------------------------------------------------
# Class for defining, parsing and evaluating the summary measures.
# Expressions sVar.exprs are evaluated in the environment of a given data.frame.
#----------------------------------------------------------------------------------
  # II) sVar naming:
  #todo 34 (sVar_evaluator, sVar.name) +0: Instead of throwing an error for un-named vector results (e.g., def.sW(A)), 
    # would be nice to extract the first is.name of the expression and use that as a name.
    # this would also eliminate the need for noname argument -> when def.sA(W1) or def.sA("W1") always set noname=TRUE 
  #todo 7 (sVar_evaluator, sVar.name) +0: Check that the resulting column names in sVar are all unique!
  # Create an option $keep.sVar.nms; When TRUE do not change the output column names for sVar mat with ncol > 1. 
  # Create an option to overwrite sVar mat colname(s) with user-provided names
  #todo 28 (sVar_evaluator) +0: Consider returning sVar.res_l instead of mat.sVar, also see if there are faster alternatives to cbind
    # (i.e., pre-allocating sVar.mat); perform benchmarks to see if there is any noticable benefit
  #todo 42 ('+..DefineSummariesClass') +0: Allow adding character vector summary measures for sVar2, s.a., 
  # def.sW(W2[[1:Kmax]]) + "netW3_sum = rowSums(W3[[1:Kmax]]"

# Useful function for testing if a name is a valid R object name:
isValidAndUnreservedName <- function(string) {
  make.names(string) == string 	# make.names converts any string into a valid R object name
}

# Wrappers for DefineSummariesClass$new(...) constructor:
#' Define Summary Measures sA and sW
#'
#' Define and store \code{tmlenet} function arguments \code{sW} and \code{sA},
#'  which are arbitrary summary measures based on treatmet \code{Anode} and baseline covariates in \code{Wnodes} 
#'  constructed from the input \code{data}.
#'  Each summary measure \code{sVar} is specified as a named R expression or a named character argument. 
#'  Separate calls to \code{def.sW/def.sA} functions can be aggregated into a single collection with '+' function, 
#'  e.g., \code{def.sW(W1=W1)+def.sW(W2=W2)}.
#'  Two special named arguments can be passed to the \code{def.sW}, \code{def.sA} functions:
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
#' @seealso \code{\link{eval.summaries}} for
#'  evaluation and validation of the summary measures,
#'  \code{\link{tmlenet}} for estimation.
#' @example tests/examples/2_defsWsA_examples.R
#' @export
def.sW <- function(...) {
  # --------------------------------------------------------------------------------------------------
  # call outside fun that parses ... and assigns empty names "" if the names attribute not set:
    # ...
    sVar.exprs <- eval(substitute(alist(...)))
    if (length(sVar.exprs)>0) { # deparse into characters when expr is.call, but keep as-is otherwise
      sVar.exprs <- lapply(sVar.exprs, function(x) if (is.character(x)) {x} else {deparse(x)})
    }
    # if not a single argument was named, names attribute of sVar.exprs will be null => add names attribute
    if (is.null(names(sVar.exprs))) names(sVar.exprs) <- rep_len("", length(sVar.exprs))
    if (length(sVar.exprs)!=0 && any(names(sVar.exprs)%in%"")) {
      message("Some summary measures were not named, automatic column name(s) will be generated during evaluation")
    }
  # end of outside fun
  # --------------------------------------------------------------------------------------------------

  sVar.exprs <- c(sVar.exprs, list(nF="nF")) # add nF node (vector with counts of friends):
  node_evaluator <- DefineSummariesClass$new(type = "sW")
  node_evaluator$set.user.env(user.env = parent.frame())
  node_evaluator$set.new.exprs(exprs_list = sVar.exprs)
  return(node_evaluator)
}

# old version:
# def.sW <- function(...) {
#   Define_sVar$new(..., type = "sW", user.env = parent.frame())
# }
# **********************************************************************
# TODO: Consider adding argument Anode to def.sA function
# **********************************************************************
#' @rdname def.sW
#' @export
def.sA <- function(...) {
  # --------------------------------------------------------------------------------------------------
  # call outside fun that parses ... and assigns empty names "" if the names attribute not set:
    # ...
    sVar.exprs <- eval(substitute(alist(...)))
    if (length(sVar.exprs)>0) { # deparse into characters when expr is.call, but keep as-is otherwise
      sVar.exprs <- lapply(sVar.exprs, function(x) if (is.character(x)) {x} else {deparse(x)})
    }
    # if not a single argument was named, names attribute of sVar.exprs will be null => add names attribute
    if (is.null(names(sVar.exprs))) names(sVar.exprs) <- rep_len("", length(sVar.exprs))
    if (length(sVar.exprs)!=0 && any(names(sVar.exprs)%in%"")) {
      message("Some summary measures were not named, automatic column name(s) will be generated during evaluation")
    }
  # end of outside fun
  # --------------------------------------------------------------------------------------------------

  node_evaluator <- DefineSummariesClass$new(type = "sA")
  node_evaluator$set.user.env(user.env = parent.frame())
  node_evaluator$set.new.exprs(exprs_list = sVar.exprs)
  return(node_evaluator)
}

is.DefineSummariesClass <- function(obj) "DefineSummariesClass" %in% class(obj)

# S3 method '+' for adding two Define_sVar objects
# Summary measure lists in both get added as c(,) into the summary measures in sVar1 object
#' @rdname def.sW
#' @param sVar1 An object returned by a call to \code{def.sW} or \code{def.sA} functions.
#' @param sVar2 An object returned by a call to \code{def.sW} or \code{def.sA} functions.
#' @export
`+.DefineSummariesClass` <- function(sVar1, sVar2) {
  assert_that(is.DefineSummariesClass(sVar1))
  assert_that(is.DefineSummariesClass(sVar2))
  assert_that(all.equal(sVar1$type, sVar2$type))
  sVar1$add.new.exprs(NewSummaries = sVar2)
  return(sVar1)
}

# ------------------------------------------------------------
# DEPRECATED PARSER (TO BE REMOVED):
# ------------------------------------------------------------
# take sVar expression index and evaluate:
parse.sVar.out <- function(sVar.idx, self, data.df) {
  sVar.expr <- self$exprs_list[[sVar.idx]]
  sVar.name <- names(self$exprs_list)[sVar.idx]
  # sVar.name <- self$sVar.expr.names[sVar.idx]
  misXreplace <- self$sVar.misXreplace[sVar.idx]
  # ******
  # eval.sVar.params <- c(self$df.names(data.df), list(misXreplace = misXreplace), list(netind_cl = self$netind_cl))
  eval.sVar.params <- c(list(self = self),
                        self$df.names(data.df), # special var "ANCHOR_ALLVARNMS_VECTOR_0" with names of already simulated vars
                        list(misXreplace = misXreplace), # replacement value for missing network covars
                        list(netind_cl = self$netind_cl),
                        list(nF = self$netind_cl$nF)
                        )
  data.env <- c(eval.sVar.params, self$node_fun, data.df)
  # ******

  if (is.character(sVar.expr)) {
    sVar.expr_call <- try(parse(text=sVar.expr)[[1]])   # parse expression into a call
    if(inherits(sVar.expr_call, "try-error")) {
      stop("error while evaluating expression: " %+% sVar.expr %+% ".\nCheck syntax specification.", call.=FALSE)
    }
  } else if (is.call(sVar.expr)){
    sVar.expr_call <- sVar.expr
    message(sVar.expr_call %+% ": sVar formula is already a parsed call")
  } else {
    stop("sVar formula class: " %+% class(sVar.expr) %+% ". Currently can't process sVar formulas that are not strings or expressions")
  }

  # print("sVar.name: "); print(sVar.name)
  # print("sVar.expr_call: "); print(sVar.expr_call)

  sVar.expr_call <- eval(substitute(substitute(e, list(Kmax = eval(self$Kmax))), list(e = sVar.expr_call))) # Replace Kmax its val
  evaled_expr <- try(eval(sVar.expr_call, envir = data.env, enclos = self$user.env)) # eval'ing expr in the envir of data.df

  # no.sVar.name <- self$sVar.noname[sVar.idx]
  no.sVar.name <- is.null(sVar.name) || (sVar.name %in% "")

  if (is.matrix(evaled_expr)) {
    if (no.sVar.name) sVar.name <- colnames(evaled_expr)
    if (!no.sVar.name && ncol(evaled_expr) > 1) sVar.name <- sVar.name %+% "." %+% (1 : ncol(evaled_expr))
    # if (no.sVar.name || ncol(evaled_expr) > 1) message(sVar.expr %+% ":
      # the result matrix is assigned the following column name(s): " %+% paste(sVar.name, collapse = ","))
  } else {
    # if (no.sVar.name) stop(sVar.expr %+% ": summary measures not using the syntax Var[[...]] can't use argument 'noname=TRUE'")
    evaled_expr <- as.matrix(evaled_expr)
  }
  colnames(evaled_expr) <- sVar.name
  return(evaled_expr)
}

# ------------------------------------------------------------------------------------------
# Standardize all names (and fill-in the empty names) according TO THE *SAME* *NAMING* *CONVENTION*;
# ....................................
# Naming conventions for evaluation results of def.sW & def.sA:
# to be completed
# ....................................
# ------------------------------------------------------------------------------------------
# *) Will modify the results of evaluation:
  # 1) If result is a vector -> make into one column matrix
  # 2) Assign the column names to each matrix using the following naming schematic:
eval.standardize.expr <- function(expr.idx, self, data.df) {

  # -------------------------------------------------------
  # First evaluate the expression result:
  # -------------------------------------------------------
  evalres <- eval.nodeform.out(expr.idx = expr.idx, self = self, data.df = data.df)
  # expression itself as string:
  expr_char <- self$exprs_list[[expr.idx]]
  # current expression name:
  expr_nm <- names(self$exprs_list)[expr.idx]
  # names of parents vars for this expression:
  expr_parents <- evalres[["par.nodes"]]

  # -------------------------------------------------------
  # no user-supplied argument name, hence need to name this expression:
  # -------------------------------------------------------
  # flag TRUE if user did not provide an argument name for self$exprs_list[expr.idx]:
  expr_noname <- (names(self$exprs_list)[expr.idx] %in% "")
  message(expr_char %+% " expr_noname? " %+% expr_noname);
  if (expr_noname && (length(expr_parents)>1)) {
    # as an alternative, can use a name paste0(expr_parents, collapse=".")
    stop("must name complex expressions that involve more than one variable: " %+% expr_char)
  } else if (expr_noname && !is.null(expr_parents) && is.character(expr_parents)) {
    message("assigning a name '" %+% expr_parents %+% "' to expression: " %+% expr_char)
    expr_nm <- expr_parents
  } else if (expr_noname) stop(expr_char%+% ": parents are null or not a character vector")

  # -------------------------------------------------------
  # convert vectors to columns, name the matrix columns according to the same naming convention:
  # -------------------------------------------------------
  # evaluation result:
  # expr_res <- evalres[["evaled_expr"]]
  # if result a vector: convert to one-col matrix assign a name: names(self$exprs_list)[expr.idx] = expr_parents
  if (is.vector(evalres[["evaled_expr"]])) {
     print(expr_char %+% ": expression result is a vector, converting to 1 col matrix and assigning a col name, " %+% expr_nm)
     expr_res <- matrix(data = evalres[["evaled_expr"]], ncol = 1)
     colnames(expr_res) <- expr_nm
     return(list(new_expr_name = expr_nm, evaled_expr = expr_res, par.nodes = evalres[["par.nodes"]]))
  # for matrix results: if column names exist (!is.null(colnames(expr_res))): DO NOTHING
  # if column names don't exist (is.null(colnames(expr_res))) or some are empty strings "":
  } else if (is.matrix(evalres[["evaled_expr"]])) {
    if (is.null(colnames(evalres[["evaled_expr"]])) || (any(colnames(evalres[["evaled_expr"]])%in%"")) || (length(expr_parents)>1)) {
      # assign names by convention: <- expr_nm%+%"."%+%c(1:ncol(expr_res))
      message(expr_char %+% ": assigning new column names to the evaluation result")
      colnames(evalres[["evaled_expr"]]) <- expr_nm%+%"."%+%c(1:ncol(evalres[["evaled_expr"]]))
    }
    return(list(new_expr_name = expr_nm, evaled_expr = evalres[["evaled_expr"]], par.nodes = evalres[["par.nodes"]]))
  } else {
    # if result is not a vector or matrix: throw an exception
    stop(expr_char%+% ": summary measure result type is " %+%class(evalres[["evaled_expr"]])%+%"; only matrix or a vector results are supported")
  }
}

## ---------------------------------------------------------------------
#' R6 class for parsing and evaluating user-specified summary measures (in \code{exprs_list})
# **********************************************************************
# TO DO: Add the description of the naming conventions for matrix summary measures
# **********************************************************************
#'
#' This R6 class can parse and evaluate (given the input data) the summary measures defined with functions 
#'  \code{\link{def.sW}} and \code{\link{def.sA}}. 
#'  The summary measure expressions (stored in \code{exprs_list}) are evaluated in the environment of the input data.frame.
#'  Note that the resulting summary measures are never stored inside this class, 
#'  the data is stored only in \code{\link{DatNet}} and \code{\link{DatNet.sWsA}} R6 classes.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{type}} - Type of the summary measure, sW or sA, determined by the calling functions \code{def.sW} or \code{def.sA}.
#' \item{\code{exprs_list}} - Deparsed list of summary measure expressions (as strings).
#' \item{\code{new_expr_names}} - Re-evaluated summary measure names, if none were provided by the user these will be 
#'  evaluated on the basis of variable names used in the expression.
#' \item{\code{sVar.names.map}} - Map of the summary measure names given by the user to actual column names in the summary measure matrix.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(type)}}{...}
#'   \item{\code{set.new.exprs(exprs_list)}}{...}
#'   \item{\code{add.new.exprs(NewSummaries)}}{...}
#'   \item{\code{eval.nodeforms(data.df, netind_cl)}}{...}
#'   \item{\code{df.names(data.df)}}{List of variables in the input data \code{data.df} gets assigned to a special 
#'   variable (\code{ANCHOR_ALLVARNMS_VECTOR_0})}
#' }
#' @importFrom assertthat assert_that
#' @export
DefineSummariesClass <- R6Class("DefineSummariesClass",
# Define_sVar_tmlenet <- R6Class("Define_sVar",
  class = TRUE,
  portable = TRUE,
  inherit = Define_sVar,
  public = list(
    type = NA,                      # "sW" or "sA" depending on which functions called the constructor: def.sW or def.sA
    exprs_list = list(),            # list of expressions in character strings, with attribute "names" set to the user-supplied names of the expressions
    new_expr_names = list(),        # re-evaluated summary measure names, if non provided by the user these will be evaluated on the basis of variable names used in the expression
    sVar.names.map = list(),        # the map between user-supplied expression (argument names) to the column names of each expression in self$exprs_list

    # exprs_evalres = list(),       # the results of evaluation of each expression in exprs_list
    # exprs_parents = list(),       # list of vectors with data.df variable names referenced in each expression, extracted by parser on each expression
    # same as above sVar.names.map
    # exprs_res_names = list(),     # list of vectors with the names of the evaluation column/vectors for each expression exprs_list[[i]] (if result is a matrix, its a vector of column names, if result is unnamed vector, its "")
    # sVar.noname = FALSE,          # (TO BE REMOVED) vector, for each TRUE sVar.expr[[idx]] ignores user-supplied name and generates names automatically
    # data.df = NULL,               # data.frame that is used for evaluation of sVar expressions (passed to get.mat.sVar)
    # Kmax = NULL,
    # netind_cl = NULL,
    # ReplMisVal0 = FALSE,          # Replace missing network VarNode values (when nF[i] < Kmax) with gvars$misXreplace (0)?
    # sVar.misXreplace = NULL,      # Replacement values for missing sVar (length(exprs_list)), either gvars$misXreplace or gvars$misval
    # sVar.noname = FALSE,          # Vector, for each TRUE sVar.expr[[idx]] ignores user-supplied name and generates names automatically
    # exprs_list = character(),     # deparsed sVar expressions (char vector)
    # sVar.expr.names = character(),# user-provided name of each sVar.expr
    # user.env = emptyenv(),        # User environment used as enclos arg to eval(sVar, enclos=)
    # mat.sVar = matrix(),        # no longer storing the sVar evaluation result

    initialize = function(type) {
      self$type <- type
      invisible(self)
    },

    # define new summary measures to be evaluated later:
    set.new.exprs = function(exprs_list) {
      # will define 1) self$exprs_list; 2) self$asis.flags; 3) self$ReplMisVal0; 4) self$sVar.misXreplace
      super$set.new.exprs(exprs_list)
      # map from self$exprs_list to variable names in each expression (column names)
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

# ---------------------------------------------------------------------------
# DEPRECATED FUNCTION (TO BE REMOVED):
# ---------------------------------------------------------------------------
    get.mat.sVar = function(data.df, netind_cl, addnFnode = NULL) {
      assert_that(is.data.frame(data.df))
      if (missing(netind_cl) && is.null(self$netind_cl)) stop("must specify netind_cl arg at least once")
      if (!missing(netind_cl)) self$netind_cl <- netind_cl
      self$Kmax <- self$netind_cl$Kmax
      # call lapply on parse.sVar.out for each sVar in sVar.expr.names -> sVar.res_l
      sVar.res_l <- lapply(seq_along(self$exprs_list), parse.sVar.out, self = self, data.df = data.df)
      names(sVar.res_l) <- names(self$exprs_list)
      if (!is.null(addnFnode)) sVar.res_l <- c(sVar.res_l, list(nF = netind_cl$mat.nF(addnFnode)))
      # SAVE THE MAP BETWEEEN EXPRESSION NAMES AND CORRESPONDING COLUMN NAMES:
      self$sVar.names.map <- lapply(sVar.res_l, colnames)
      names(self$sVar.names.map) <- names(self$exprs_list)
      mat.sVar <- do.call("cbind", sVar.res_l)
      return(mat.sVar)
    },

   eval.nodeforms = function(data.df, netind_cl) {
      assert_that(is.data.frame(data.df))
      if (missing(netind_cl) && is.null(self$netind_cl)) stop("must specify netind_cl arg at least once")
      if (!missing(netind_cl)) self$netind_cl <- netind_cl
      self$Kmax <- self$netind_cl$Kmax
      self$Nsamp <- nrow(data.df)

      # EVALUATE THE EXPRESSION AND THEN STANDARDIZE ALL NAMES ACCORDING TO ONE NAMING CONVENTION:
      sVar.res_l <- lapply(seq_along(self$exprs_list), eval.standardize.expr, self = self, data.df = data.df)
      # sVar.res_l <- lapply(seq_along(self$exprs_list), eval.nodeform.out, self = self, data.df = data.df)
      self$new_expr_names <- lapply(sVar.res_l, '[[', 'new_expr_name')
      names(self$new_expr_names) <- unlist(self$new_expr_names)

      print("old names(self$exprs_list): "); print(names(self$exprs_list))
      print("new names self$new_expr_names: "); print(names(self$new_expr_names))

      names(sVar.res_l) <- names(self$new_expr_names)
      names(self$exprs_list) <- names(self$new_expr_names)

      # assign self$sVar.names.map based on newly standardized summary names:
      self$sVar.names.map <- lapply(sVar.res_l, function(x) colnames(x[["evaled_expr"]]))
      names(self$sVar.names.map) <- names(self$new_expr_names)
      print("self$sVar.names.map: "); print(str(self$sVar.names.map))

      mat.sVar <- do.call("cbind", lapply(sVar.res_l, function(x) x[["evaled_expr"]]))
      print("eval result mat after standardizing: "); print(head(mat.sVar))
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