### --- Test setup ---
 
if(FALSE) {
  library("RUnit")
  library("roxygen2")
  library("devtools")
  # setwd(".."); setwd(".."); getwd()
  document()
  setwd("..")
  install("tmlenet", build_vignettes=FALSE) # INSTALL W/ devtools:

  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf tmlenet")  # just create the pdf manual from help files

  # CHECK AND BUILD PACKAGE:
  # setwd(".."); setwd(".."); getwd()
  # devtools::check() # runs check with devtools

  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  devtools::build(args = "--compact-vignettes") # build package tarball compacting vignettes
  # devtools::build() # build package tarball
  setwd("..")
  system("R CMD check --as-cran tmlenet_0.0.9.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes tmlenet") # check without building the pdf manual and not building vignettes
      ## system("R CMD build tmlenet --no-build-vignettes")
      ## system("R CMD build tmlenet")  
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./tmlenet_0.0.9.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(tmlenet)
  # tmlenet:::addvectorfcn("poisson")
  # tmlenet:::debug_set() # SET TO DEBUG MODE
  # tmlenet:::debug_off() # SET DEBUG MODE OFF

  # TEST COVERATE:
  # if your working directory is in the packages base directory
  # package_coverage()
  # or a package in another directory
  # cov <- package_coverage("tmlenet")
  # view results as a data.frame
  # as.data.frame(cov)
  # zero_coverage() can be used to filter only uncovered lines.
  # zero_coverage(cov)
}

psi_RDs_DAG2a <- NULL
psi_RDs_DAG2b <- NULL

sample_checks <- function() {   # doesn't run, this is just to show what test functions can be used
  print("Starting tests...")
    checkTrue(1 < 2, "check1")     ## passes fine
     ## checkTrue(1 > 2, "check2")  ## appears as failure in the test protocol
     v <- 1:3
     w <- 1:3
     checkEquals(v, w)               ## passes fine
     names(v) <- c("A", "B", "C")
     ## checkEquals(v, w)            ## fails because v and w have different names
     checkEqualsNumeric(v, w)        ## passes fine because names are ignored
     x <- rep(1:12, 2)
     y <- rep(0:1, 12)
     res <- list(a=1:3, b=letters, LM=lm(y ~ x))
     res2 <- list(a=seq(1,3,by=1), b=letters, LM=lm(y ~ x))
     checkEquals( res, res2)        ## passes fine
     checkIdentical( res, res)
     checkIdentical( res2, res2)
     ## checkIdentical( res, res2)  ## fails because element 'a' differs in type
     fun <- function(x) {
       if(x)
       {
        stop("stop conditions signaled")
       }
       return()
     }
     checkException(fun(TRUE))      ## passes fine
     ## checkException(fun(FALSE))  ## failure, because fun raises no error
     checkException(fun(TRUE), silent=TRUE)
     ##  special constants
     ##  same behaviour as for underlying base functions
     checkEquals(NA, NA)
     checkEquals(NaN, NaN)
     checkEquals(Inf, Inf)
     checkIdentical(NA, NA)
     checkIdentical(NaN, NaN)
     checkIdentical(-Inf, -Inf)
}

`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))

test.bugfixes <- function() {
}

test.plotting <- function() {
}



test.node <- function() {
}


test.faster_tolongdata <- function() {
    #-------------------------------------------------------------
    # Testing the converter to long format that is based on package data.table
    #-------------------------------------------------------------
    # t_end <- 16
    # library(simcausal)

    # D <- DAG.empty()
    # D <- D + node("L2", t=0,        distr="rbern", prob=0.05, order=1)
    # D <- D + node("L1", t=0,        distr="rbern", prob=ifelse(L2[0]==1,0.5,0.1), order=2)
    # D <- D + node("A1", t=0,        distr="rbern", prob=ifelse(L1[0]==1 & L2[0]==0, 0.5, ifelse(L1[0]==0 & L2[0]==0, 0.1, ifelse(L1[0]==1 & L2[0]==1, 0.9, 0.5))), order=3)
    # D <- D + node("A2", t=0,        distr="rbern", prob=0, order=4, EFU=TRUE)
    # D <- D + node("Y",  t=0,        distr="rbern", prob=plogis(-6.5 + L1[0] + 4*L2[0] + 0.05*I(L2[0]==0)), order=5, EFU=TRUE)
    # D <- D + node("L2", t=1:t_end,  distr="rbern", prob=ifelse(A1[t-1]==1, 0.1, ifelse(L2[t-1]==1, 0.9, min(1,0.1 + t/16))), order=6+4*(0:(t_end-1)))
    # D <- D + node("A1", t=1:t_end,  distr="rbern", prob=ifelse(A1[t-1]==1, 1, ifelse(L1[0]==1 & L2[0]==0, 0.3, ifelse(L1[0]==0 & L2[0]==0, 0.1, ifelse(L1[0]==1 & L2[0]==1, 0.7, 0.5)))), order=7+4*(0:(t_end-1)))
    # D <- D + node("A2", t=1:t_end,  distr="rbern", prob=plogis(-3.5 + 0.5*A1[t]+0.5*L2[t]), order=8+4*(0:(t_end-1)), EFU=TRUE) # informative censoring

    # # this takes longer (6 sec longer for 1Mil obs)
    # # D <- D + node("Y",  t=1:t_end,  distr="rbern", prob=plogis(-6.5 + L1[0] + 4*L2[t] + 0.05*sum(I(L2[0:t]==rep(0,(t+1))))), order=9+4*(0:(t_end-1)), EFU=TRUE)
    # D <- D + node("Y",  t=1:t_end,  distr="rbern", prob=plogis(-6.5 + L1[0] + 4*L2[t] + 0.05*(sum(L2[0:t])==0)), order=9+4*(0:(t_end-1)), EFU=TRUE)
    # lDAG3 <- set.DAG(D)

    # #-------------------------------------------------------------
    # # Adding dynamic actions (indexed by a real-valued parameter)
    # #-------------------------------------------------------------
    # # act_t0_theta <- node("A1",t=0, distr="rbern", prob=ifelse(L2[0] >= theta,1,0))
    # act_t0_theta <- node("A1",t=0, distr="rbern", prob=ifelse(L2[0] >= theta,1,0))
    # act_tp_theta <- node("A1",t=1:t_end, distr="rbern", prob=ifelse(A1[t-1]==1,1,ifelse(L2[t] >= theta,1,0)))
    # act_NoCens <- node("A2",t=0:t_end, distr="rbern", prob=0)
    # actionnodes <- c(act_t0_theta, act_tp_theta, act_NoCens)
    # D <- lDAG3 + action("A1_th0", nodes=actionnodes, theta=0)
    # D <- D + action("A1_th1", nodes=actionnodes, theta=1)

    # #-------------------------------------------------------------
    # # Testing conversion of observed data to long format 
    # #-------------------------------------------------------------
    # # NO carry forward imputation:
    # O_dat_df <- simobs(D, n=500, rndseed = 123)
    # system.time(O_dat_long <- DF.to.long(O_dat_df))
    # system.time(O_dat_long_DT <- DF.to.longDT(O_dat_df))

    # checkIdentical(O_dat_long$ID, O_dat_long_DT$ID)
    # checkIdentical(O_dat_long$L1_0, O_dat_long_DT$L1_0)
    # checkIdentical(O_dat_long$t, O_dat_long_DT$t)
    # checkIdentical(O_dat_long$L2, O_dat_long_DT$L2)
    # checkIdentical(O_dat_long$A1, O_dat_long_DT$A1)
    # checkIdentical(O_dat_long$A2, O_dat_long_DT$A2)
    # checkIdentical(O_dat_long$Y, O_dat_long_DT$Y)

    # # With carry forward imputation of Y (vs 1):
    # O_dat_df <- simobs(D, n=500, rndseed = 123)
    # O_dat_LTCF <- doLTCF(data=O_dat_df, LTCF="Y")
    # system.time(O_dat_long_LTCF_v1 <- DF.to.long(O_dat_LTCF))
    # system.time(O_dat_long_DT_LTCF_v1 <- DF.to.longDT(O_dat_LTCF))

    # checkIdentical(O_dat_long_LTCF_v1$ID, O_dat_long_DT_LTCF_v1$ID)
    # checkIdentical(O_dat_long_LTCF_v1$L1_0, O_dat_long_DT_LTCF_v1$L1_0)
    # checkIdentical(O_dat_long_LTCF_v1$t, O_dat_long_DT_LTCF_v1$t)
    # checkIdentical(O_dat_long_LTCF_v1$L2, O_dat_long_DT_LTCF_v1$L2)
    # checkIdentical(O_dat_long_LTCF_v1$A1, O_dat_long_DT_LTCF_v1$A1)
    # checkIdentical(O_dat_long_LTCF_v1$A2, O_dat_long_DT_LTCF_v1$A2)
    # checkIdentical(O_dat_long_LTCF_v1$Y, O_dat_long_DT_LTCF_v1$Y)

    # # With carry forward imputation of Y (vs 2):
    # O_dat_df_LTCF <- simobs(D, n=500, LTCF="Y", rndseed = 123)
    # system.time(O_dat_long_LTCF <- DF.to.long(O_dat_df_LTCF))
    # system.time(O_dat_long_DT_LTCF <- DF.to.longDT(O_dat_df_LTCF))

    # checkIdentical(O_dat_long_LTCF$ID, O_dat_long_DT_LTCF$ID)
    # checkIdentical(O_dat_long_LTCF$L1_0, O_dat_long_DT_LTCF$L1_0)
    # checkIdentical(O_dat_long_LTCF$t, O_dat_long_DT_LTCF$t)
    # checkIdentical(O_dat_long_LTCF$L2, O_dat_long_DT_LTCF$L2)
    # checkIdentical(O_dat_long_LTCF$A1, O_dat_long_DT_LTCF$A1)
    # checkIdentical(O_dat_long_LTCF$A2, O_dat_long_DT_LTCF$A2)
    # checkIdentical(O_dat_long_LTCF$Y, O_dat_long_DT_LTCF$Y)

    # #-------------------------------------------------------------
    # # Testing conversion of full data to long format (with carry forward imputation)
    # #-------------------------------------------------------------
    # X_dat <- simfull(A(D), n=500, rndseed = 123)
    # attributes(X_dat[[1]])$node_nms

    # system.time(X_dat_l <- lapply(X_dat, DF.to.long))
    # system.time(X_dat_lDT <- lapply(X_dat, DF.to.longDT))

    # checkIdentical(X_dat_l[["A1_th0"]]$ID, X_dat_lDT[["A1_th0"]]$ID)
    # checkIdentical(X_dat_l[["A1_th0"]]$L1_0, X_dat_lDT[["A1_th0"]]$L1_0)
    # checkIdentical(X_dat_l[["A1_th0"]]$t, X_dat_lDT[["A1_th0"]]$t)
    # checkIdentical(X_dat_l[["A1_th0"]]$L2, X_dat_lDT[["A1_th0"]]$L2)
    # checkIdentical(X_dat_l[["A1_th0"]]$A1, X_dat_lDT[["A1_th0"]]$A1)
    # checkIdentical(X_dat_l[["A1_th0"]]$A2, X_dat_lDT[["A1_th0"]]$A2)
    # checkIdentical(X_dat_l[["A1_th0"]]$Y, X_dat_lDT[["A1_th0"]]$Y)

    # checkIdentical(X_dat_l[["A1_th1"]]$ID, X_dat_lDT[["A1_th1"]]$ID)
    # checkIdentical(X_dat_l[["A1_th1"]]$L1_0, X_dat_lDT[["A1_th1"]]$L1_0)
    # checkIdentical(X_dat_l[["A1_th1"]]$t, X_dat_lDT[["A1_th1"]]$t)
    # checkIdentical(X_dat_l[["A1_th1"]]$L2, X_dat_lDT[["A1_th1"]]$L2)
    # checkIdentical(X_dat_l[["A1_th1"]]$A1, X_dat_lDT[["A1_th1"]]$A1)
    # checkIdentical(X_dat_l[["A1_th1"]]$A2, X_dat_lDT[["A1_th1"]]$A2)
    # checkIdentical(X_dat_l[["A1_th1"]]$Y, X_dat_lDT[["A1_th1"]]$Y)

    # BENCHMARKING for 50K: gain of x5.4 factor
    # old convert to long
    # user  system elapsed 
    # 13.943   1.783  15.766 
    # new DF.to.longDT
    # user  system elapsed 
    # 2.564   0.370   2.935 

    # BENCHMARKING for 500K: gain of x5.4 factor
    # old convert to long
    # user  system elapsed 
    # 140.378  18.398 158.853     
    # new DF.to.longDT    
    # user  system elapsed 
    # 28.753   4.092  32.844 

    # CONVERTING BACK TO WIDE FORMAT:
    # ## convert long to wide using dcast from data.table
    # dcast.data.table(data, formula, fun.aggregate = NULL, 
    #     ..., margins = NULL, subset = NULL, fill = NULL, 
    #     drop = TRUE, value.var = guess(data),
    #     verbose = getOption("datatable.verbose"))


    #-------------------------------------------------------------
    # BENCHMARKING current package rowSums with data.table version - abandoned for now
    #-------------------------------------------------------------
    # t_pts <- 0:16
    # t_vec <- "_"%+%(t_pts)
    # (L2_names <- "L2"%+%t_vec)
    # # library(data.table)
    # system.time( X_dat_1 <- simfull(A(D)[1], n=1000000, LTCF="Y", rndseed = 123)[[1]])
    # X_dat_1 <- X_dat_1[,-1]
    # nrow(X_dat_1)
    # colnames(X_dat_1)
    # head(X_dat_1)
    # X_dat_1_DT = data.table(X_dat_1)
    # setkey(X_dat_1_DT,ID)
    # L2_idx <- which(names(X_dat_1_DT)%in%L2_names)
    # ID_idx <- which(names(X_dat_1_DT)%in%"ID")
    # ncol(X_dat_1_DT)
    # # COMPARING data.table to current data.rame rowSums
    # # fast version 1  of row sums 1
    # system.time(
    # X_dat_1_DT[, SumRow := rowSums(.SD), .SDcols = L2_idx]
    # )
    # user  system elapsed 
    # 0.139   0.030   0.181     

    # # faster version 2 of row sums
    # system.time(
    # X_dat_1_DT[, SumRow2 := Reduce(`+`, .SD), .SDcol = L2_idx]
    # )
    #  user  system elapsed 
    # 0.075   0.027   0.120

    # # version 3 using set
    # # i  In set(), integer row numbers to be assigned value.
    # # NULL represents all rows more efficiently than creating a vector such as 1:nrow(x).
    # # j In set(), integer column number to be assigned value.    
    # # x - DT, i - row indx, j - col indx, value - new val
    # set(x, i=NULL, j, value)
    # system.time(for (i in 1:nrow(X_dat_1_DT)) set(X_dat_1_DT,as.integer(i),"SumRow3",i))

    # # current version of row sums (slowest)
    # system.time(newVar_vold <- rowSums(X_dat_1[,L2_names]))
    # # COMPARING data.table to current data.rame rowSums with some column operations
    # system.time(X_dat_1_DT[, SumRow := rowSums(.SD==c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), .SDcols = L2_idx])
    # system.time(newVar_vold <- rowSums(I(X_dat_1[,L2_names] == cbind(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))))

}

