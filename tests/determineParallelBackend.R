###################################################################################
### setParallelBackend.R
### Script runs through logic to determine whether doParallel (forking) or doMPI (cluster) should
### be used to register a backend.
### Will fail if either doParallel or doMPI are not installed.
### After finishing have to terminate the MPI cluster with:
### if (useMPI) {
###   cat("closing cluster and finalizing MPI...\n")
###   closeCluster(cl)
###   mpi.finalize()
### }
###################################################################################
nHosts <- as.numeric(Sys.getenv('NHOSTS'))
nCores <- as.numeric(Sys.getenv('NSLOTS'))
if (is.na(nHosts) || nHosts == 1) { # running locally or on a single host
  useMPI <- FALSE
  require(doParallel)
  # check to see if running locally, so you can grab all cores
  if (is.na(nCores)) nCores <- detectCores() 
  registerDoParallel(cores=nCores)
  cat("running script inside a single node (", nCores, "CPUs) via doParallel...\n")
} else {
  useMPI <- TRUE
  cat("running script across", nHosts, "nodes with", nCores, "cores via MPI...\n")
  require(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)
  exportDoMPI(cl,varlist=ls())
}
# # testing parallel execution:
# test.results <- foreach(id=icount(100)) %dopar%
# {
#   res <- rnorm(10^7)
#   res <- rnorm(10^7)
#   return(c(mean(res), id))
# }