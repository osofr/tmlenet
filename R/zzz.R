# .onAttach <- function(...) {
# 	packageStartupMessage('simcausal')
# 	packageStartupMessage('Version: ', utils::packageDescription('simcausal')$Version)
# 	packageStartupMessage('Package created on ', utils::packageDescription('simcausal')$Date, '\n')
# 	packageStartupMessage('Please note this package is still in its early stages of development. Check for updates and report bugs at http://github.com/osofr/simcausal.', '\n')
# 	packageStartupMessage('To see the vignette use vignette("simcausal_vignette", package="simcausal"). To see all available package documentation use help(package = "simcausal") and ?simcausal.', '\n')
# 	packageStartupMessage('To see the latest updates for this version, use news(package = "simcausal").', '\n')
# }

.onAttach <- function(...) {
  packageStartupMessage("tmlenet")
  packageStartupMessage("The tmlenet package is still in beta testing. Interpret results with caution.")
}


`%+%` <- function(a, b) paste0(a, b)

gvars <- new.env(parent = emptyenv())
gvars$verbose <- FALSE      # verbose mode (print all messages)
gvars$misval <- NA_integer_ # the default missing value for observations
# gvars$misval <- 1L
# gvars$misval <- -.Machine$integer.max
gvars$misXreplace <- 0L     # the default replacement value for misval that appear in the design matrix
gvars$tolerr <- 10^-12      # tolerance error: assume for abs(a-b) < gvars$tolerr => a = b
gvars$maxncats <- 5L        # max number of unique categories a categorical variable sA[j] can have. If sA[j] has more it is automatically considered continuous
# gvars$maxncats <- 10L
gvars$nbins <- 15L          # default n bins for continous var
# gvars$nbins <- 50L        # default n bins for continous var
gvars$max_nperbin <- 500L   # max number of observations per 1 bin (for data-adaptive bin defns) - CONSIDER REVERSING THIS BACK TO nbins alone
gvars$sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")

# returns a function (alternatively a call) that tests for missing values in sA / sW
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
set.nbins <- function(nbins) {
  old.nbins <- gvars$nbins
  gvars$nbins <- nbins
  invisible(old.nbins)
}
set.maxncats <- function(maxncats) {
  old.maxncats <- gvars$maxncats
  gvars$maxncats <- maxncats
  invisible(old.maxncats)
}
gvars$misfun <- testmisfun()