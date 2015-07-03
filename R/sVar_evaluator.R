
#----------------------------------------------------------------------------------
# Class for defining, parsing and evaluating the summary measures.
# Expressions sVar.exprs are evaluated in the environment of a given data.frame.
#----------------------------------------------------------------------------------

# TO FINISH SUMMARY MEASUREs parser need to:
  # I) sVar output class (d.f. or mat):
    # The only reason for d.f. is for subset eval., which can be re-written to work with matrices...
    # Alternatively, consider using dplyr::as_data_frame() or dplyr::data_frame() for faster conversion to d.f.s
    # Will require converting each ncol>1 matrix result of sVar[j] into a list of separate columns...

  # II) sVar naming:
  # #todo 34 (sVar_evaluator, sVar.name) +0: Instead of throwing an error for un-named vector results (i.e., def.sW.g0(A)), would be nice to extract the first is.name of the expression and use that as a name.
  # #todo 7 (sVar_evaluator, sVar.name) +0: Check that the resulting column names in sVar are all unique!
  # #todo 3 (sVar_evaluator, sVar.name) +0: Allow overwriting names for some ncol(expr_res) > 1 in sVar.res_l. Allow using naturally assigned names for some ncol(expr_res) > 1 in sVar.res_l
  # Create an option $keep.sVar.nms; When TRUE do not change the output column names for sVar mat with ncol > 1. Create an option to overwrite sVar mat colname(s) with user-provided names
  # III) Other:
  # #todo 32 (sVar_evaluator, UI) +0: Need a diagnostic tool that will evaluate and return the result of the summary measures applied to user Odata data.frame...
  # #todo 6 (sVar_evaluator, NetInd_k) +1: TEST NetInd_k, generating sA under g_star sampled data.df: nrow(NetInd_k) = n, while nrow(data.df) = n*p  => COULD BE A PROBLEM. NEED TO CHECK WORKS AS INTENDED.

# Useful function for testing if a name is a valid R object name:
isValidAndUnreservedName <- function(string) {
  make.names(string) == string 	# make.names converts any string into a valid R object name
}

# Wrappers for DefineEval.sVar$new(...) constructor:
def.sW <- function(..., type.g0 = TRUE) if (type.g0) {def.sW.g0(...)} else {def.sW.gstar(...)}
def.sA <- function(...) DefineEval.sVar$new(..., type = "sA", user.env = parent.frame())
def.sW.g0 <- function(...) DefineEval.sVar$new(..., type = "sW.g0", user.env = parent.frame())
def.sW.gstar <- function(...) DefineEval.sVar$new(..., type = "sW.gstar", user.env = parent.frame())

# take sVar expression index and evaluate:
parse.sVar.out <- function(sVar.idx, self) {
  sVar.expr <- self$sVar.exprs[[sVar.idx]]
  sVar.name <- self$sVar.expr.names[sVar.idx]
  # ******
  eval.sVar.params <- c(self$df.names(self$data.df), list(misXreplace = self$misXreplace), list(netind_cl = self$netind_cl))
  data.env <- c(eval.sVar.params, self$node_fun, self$data.df)
  # ******
  if (is.character(sVar.expr)) {
    sVar.expr_call <- try(parse(text=sVar.expr)[[1]]) 	# parse expression into a call
    if(inherits(sVar.expr_call, "try-error")) {
      stop("error while evaluating expression: " %+% sVar.expr %+% ".\nCheck syntax specification.", call.=FALSE)
    }
  } else if (is.call(sVar.expr)){
    sVar.expr_call <- sVar.expr
    warning(sVar.expr_call %+% ": sVar formula is already a parsed call")
  } else {
    stop("sVar formula class: " %+% class(sVar.expr) %+% ". Currently can't process sVar formulas that are not strings or calls.")
  }
  sVar.expr_call <- eval(substitute(substitute(e, list(Kmax = eval(self$Kmax))), list(e = sVar.expr_call))) # Replace Kmax its val
  evaled_expr <- try(eval(sVar.expr_call, envir = data.env, enclos = self$user.env)) # eval'ing expr in the envir of data.df

  no.sVar.name <- is.null(sVar.name) || (sVar.name %in% "")

  if (is.matrix(evaled_expr)) {
    if (no.sVar.name) sVar.name <- colnames(evaled_expr)
    if (!no.sVar.name && ncol(evaled_expr) > 1) sVar.name <- sVar.name %+% "." %+% (1 : ncol(evaled_expr))
    if (no.sVar.name || ncol(evaled_expr) > 1) message(sVar.expr %+% ": the result matrix is assigned the following column name(s): " %+% paste(sVar.name, collapse = ","))
  } else {
    if (no.sVar.name) stop(sVar.expr %+% ": summary measures not defined with Var[[...]] must be named.")
    evaled_expr <- as.matrix(evaled_expr)
  }
  colnames(evaled_expr) <- sVar.name
  return(evaled_expr)
}


DefineEval.sVar <- R6Class("DefineEval.sVar",
  class = TRUE,
  portable = TRUE,
  public = list(
    user.env = emptyenv(),        # user environment to be used as enclos arg to eval(sVar)
    data.df = NULL,               # data.frame that is used for evaluation of sVar expressions
    misXreplace = NULL,           # = gvars$misXreplace by default
    netind_cl = NULL,
    Kmax = NULL,
    type = NULL,                  # sW.g0, sW.gstar or sA
    sVar.exprs = character(),     # deparsed sVar expressions (char vector)
    sVar.expr.names = character(),# user-provided name of each sVar.expr
    mat.sVar = matrix(),        # matrix storing evaluated summary measures sVar

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
            # opting for replace on entire netVars_eval, do benchmarks later
            # netVars_eval[is.na(netVars_eval[, neti]), neti] <- env$misXreplace
          }
        }
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
      message("Defined the following summary expressions: "); print(self$sVar.exprs)
      self$sVar.expr.names <- names(self$sVar.exprs)
      if (length(self$sVar.exprs) != 0 && (is.null(self$sVar.expr.names) || any(self$sVar.expr.names==""))) {
        # stop("please specify a name for each sVar expression")
        message("Some summary measures were not named, automatic column name(s) will be generated during evaluation")
      }
      class(self) <- c(class(self), type)	# for later S3 dispatch
      invisible(self)
    },

    get.mat.sVar = function(data.df, netind_cl, misXreplace, addnFnode = NULL) {
      assert_that(is.data.frame(data.df))
      self$data.df <- data.df

      if (missing(netind_cl) && is.null(self$netind_cl)) stop("must specify netind_cl arg at least once")
      if (!missing(netind_cl)) self$netind_cl <- netind_cl
      self$Kmax <- self$netind_cl$Kmax

      if (!missing(misXreplace)) {
        self$misXreplace <- misXreplace
      } else {
        self$misXreplace <- gvars$misXreplace
      }

      # call lapply on parse.sVar.out for each sVar in sVar.expr.names -> sVar.res_l
      sVar.res_l <- lapply(seq_along(self$sVar.exprs), parse.sVar.out, self = self)
      names(sVar.res_l) <- self$sVar.expr.names
      if (!is.null(addnFnode)) sVar.res_l <- c(sVar.res_l, list(nF = netind_cl$mat.nF(addnFnode)))
      # print("sVar.res_l: "); print(sVar.res_l)
      # print("data.frame(sVar.res_l): "); print(data.frame(sVar.res_l))
      # print("dplyr::as_data_frame(sVar.res_l): "); print(dplyr::as_data_frame(sVar.res_l))

      # #todo 28 (sVar_evaluator) +0: Consider returning sVar.res_l instead of mat.sVar, also see if there are faster alternatives to cbind (i.e., pre-allocating sVar.mat); perform benchmarks to see if there is any noticable benefit
      mat.sVar <- do.call("cbind", sVar.res_l)
      return(mat.sVar)
      # self$mat.sVar <- do.call("cbind", sVar.res_l)
      # invisible(self)      
    },

    df.names = function(data.df) { # list of variable names from data.df with special var name (ANCHOR_ALLVARNMS_VECTOR_0)
      return(list(ANCHOR_ALLVARNMS_VECTOR_0 = colnames(data.df)))
      # allvarnms <- list(ANCHOR_ALLVARNMS_VECTOR_0 = vector())
      # allvarnms[["ANCHOR_ALLVARNMS_VECTOR_0"]] <- append(allvarnms[["ANCHOR_ALLVARNMS_VECTOR_0"]], colnames(data.df))
      # allvarnms[["ANCHOR_ALLVARNMS_VECTOR_0"]] <- unique(allvarnms[["ANCHOR_ALLVARNMS_VECTOR_0"]])
      # return(allvarnms)
    }
  ),

  active = list(
    placeholder = function() {}
  ),

  private = list(
    privplaceholder = function() {}
  )
)