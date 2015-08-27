
#-----------------------------------------------------------------------------
# Global State Vars (can be controlled globally with options(tmlenet.optname = ))
#-----------------------------------------------------------------------------
gvars <- new.env(parent = emptyenv())
gvars$verbose <- FALSE      # verbose mode (print all messages)
gvars$opts <- list()        # named list of package options that is controllable by the user (tmlenet_options())
gvars$misval <- NA_integer_ # the default missing value for observations (# gvars$misval <- -.Machine$integer.max)
gvars$misXreplace <- 0L     # the default replacement value for misval that appear in the design matrix
gvars$tolerr <- 10^-12      # tolerance error: assume for abs(a-b) < gvars$tolerr => a = b
gvars$sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")

setopt <- function(optname, val) {
  opt <- gvars$opts
  if (!(optname %in% (names(opt)))) stop(optname %+% ": this options does not exist")
  old.optval <- opt[[optname]]
  opt[[optname]] <- val
  invisible(old.optval)
}

getopt <- function(optname) {
  opt <- gvars$opts
  if (!(optname %in% (names(opt)))) stop(optname %+% ": this options does not exist")
  return(opt[[optname]])
}

#' Print Current Option Settings for \code{tmlenet}
#' @return Invisibly returns a list of \code{tmlenet} options.
#' @seealso \code{\link{tmlenet_options}}
#' @export
print_tmlenet_opts <- function() {
  print(gvars$opts)
  invisible(gvars$opts)
}

#' Setting Options for \code{tmlenet}
#'
#' Additional options that control the estimation algorithm in \code{tmlenet} package
#' @param useglm Set to \code{FALSE} to estimate with \code{speedglm::speedglm.wfit} and \code{TRUE} for \code{glm::glm.fit}.
#' @param nbins Set the default number of bins to use when binning a continous outcome variable.
#' @param maxncats Max number of unique categories a categorical variable \code{sA[j]} can have. 
#' If \code{sA[j]} has more it is automatically considered continuous.
#' @param binByMass Define bin cutoffs for a continuous outcome using even mass distribution? 
#' Setting to \code{TRUE} allows for data-adaptive bin cutoff definitions, 
#' where each bin will be defined so that it contains an approximately the same number of observations across all bins.
#' @param poolContinVar Set to \code{TRUE} for fitting a pooled regression which pools bin indicators across all bins.
#' When fitting a model for binirized continuous outcome, set to \code{TRUE} 
#' for pooling bin indicators across several bins into one outcome regression?
#' @param maxNperBin Max number of observations per 1 bin for a continuous outcome (only applies when \code{binByMass=TRUE})
#' @return Invisibly returns a list with old option settings.
#' @seealso \code{\link{print_tmlenet_opts}}
#' @export
# alpha = 0.05,
# gbound = 0.005, 
# family = "binomial", # NOT IMPLEMENTED YET
# n_MCsims = ceiling(sqrt(nrow(data))),
# onlyTMLE_B = TRUE,
# f_g0 = NULL
tmlenet_options <- function(useglm = FALSE, nbins = 15, maxncats = 5, binByMass = TRUE, poolContinVar = FALSE, maxNperBin = 1000) {
  # nbins = 50L, # maxncats = 10L
  old.opts <- gvars$opts
  opts <- list(
    useglm = useglm,
    nbins = nbins,
    maxncats = maxncats,
    binByMass = binByMass,
    poolContinVar = poolContinVar,
    maxNperBin = maxNperBin
  )
  gvars$opts <- opts
  invisible(old.opts)
}

# returns a function (alternatively a call) that tests for missing values in (sA,sW)
testmisfun <- function() {
  if (is.na(gvars$misval)) {
    return(is.na)
  } else if (is.null(gvars$misval)){
    return(is.null)
  } else if (is.integer(gvars$misval)) {
    return(function(x) {x==gvars$misval})
  } else {
    return(function(x) {x%in%gvars$misval})
  }
}

get.misval <- function() {
  gvars$misfun <- testmisfun()
  gvars$misval
}

set.misval <- function(gvars, newmisval) {
  oldmisval <- gvars$misval
  gvars$misval <- newmisval
  gvars$misfun <- testmisfun()    # EVERYTIME gvars$misval HAS CHANGED THIS NEEDS TO BE RESET/RERUN.
  invisible(oldmisval)
}
gvars$misfun <- testmisfun()

# Allows tmlenet functions to use e.g., getOption("tmlenet.verbose") to get verbose printing status
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.tmlenet <- list(
    tmlenet.verbose = gvars$verbose
  )

  tmlenet_options()

  toset <- !(names(op.tmlenet) %in% names(op))
  if(any(toset)) options(op.tmlenet[toset])

  invisible()
}

.onAttach <- function(...) {
  packageStartupMessage('tmlenet')
  packageStartupMessage('The tmlenet package is still in beta testing. Interpret results with caution.')
  #   packageStartupMessage('Version: ', utils::packageDescription('tmlenet')$Version)
  #   packageStartupMessage('Package created on ', utils::packageDescription('tmlenet')$Date, '\n')
  #   packageStartupMessage('Please note this package is still in its early stages of development.
   # Check for updates and report bugs at http://github.com/osofr/tmlenet.', '\n')
  #   packageStartupMessage('To see the vignette use vignette("tmlenet_vignette", package="tmlenet"). 
  # To see all available package documentation use help(package = "tmlenet") and ?tmlenet.', '\n')
  #   packageStartupMessage('To see the latest updates for this version, use news(package = "tmlenet").', '\n')
}











