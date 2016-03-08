`%+%` <- function(a, b) paste0(a, b)
# ---------------------------------------------------------------------------------
# TEST SET 1
# SIMPLE NETWORK TMLE
# ---------------------------------------------------------------------------------

test.tmleneterrrors <- function() {
  # Set x% of community to A=1 (returns probability P(A=1))
  f.A_x <- function(data, x, ...) rbinom(n = nrow(data), size = 1, prob = x[1])
  # Deterministically set every A=0
  f.A_0 <- function(data, ...) f.A_x(data, 0, ...)
  # Deterministically set every A=1
  f.A_1 <- function(data, ...) f.A_x(data, 1, ...)

  # incorrect intervention f_gstar1:
  f.A_wrong <- function(x, ...) 1

  #***************************************************************************************
  # EXAMPLE WITH SIMULATED DATA FOR 6 FRIENDS AND 3 W's (OLD SIMULATION 3)
  #***************************************************************************************
  data(df_netKmax6) # Load the network data
  Kmax <- 6 # Max number of friends in the network

  #----------------------------------------------------------------------------------
  # Example 1. Mean population outcome under deterministic intervention A=0 with 6 friends
  # Intercept based TMLE
  #----------------------------------------------------------------------------------
  options(tmlenet.verbose = FALSE)

  sW <- def_sW(netW2 = W2[[1:Kmax]]) +
        def_sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

  sA <- def_sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
        def_sA(netA = A[[0:Kmax]])

  # Test for non-existing predictors in hform.g0/hform.gstar:
  checkException(
    res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + netW3_sum + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    Anodes = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for non-existing predictors in hform.g0/hform.gstar:
  checkException(
    res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anodes = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for non-existing outcomes in hform.g0/hform.gstar:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "sum.netW3 ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netW3 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anodes = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for different outcomes in hform.g0/hform.gstar (non-existing for hform.gstar):
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netW3 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anodes = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for existing but different outcomes in hform.g0/hform.gstar:
  checkException(
       res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netAW2 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anodes = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Throws exception since netW3_sum, sum_1mAW2_nets from Qform don't exist:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Qform = " blah ~ netW3_sum + sum_1mAW2_nets",
                  hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",

                  Anodes = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Throws exception since Ynode arg is omitted, but blah from Qform LHS doesn't exist:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Qform = " blah ~ sum.netW3 + sum.netAW2",
                  hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",
                  Anodes = "A",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = f.A_0, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1)))

  # Throws exception when f_gstar1 function doesn't have argument "data":
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Anodes = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = f.A_wrong, sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1)))
  # Throws an exception when f_gstar1 is a vector of 1<length<n:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Anodes = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = rep(1L,50), sW = sW, sA = sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1)))
}


test.tmlenetf.gen.A.star <- function() {
  require('assertthat')
  Anodes <- "A"%+%1
  test.dat <- matrix(rnorm(50), nrow=10, ncol =5)
  abar0 <- 0
  abar1 <- rep(0,10)
  abar2 <- function(data, ...) 1
  abar3 <- function(data, ...) rep(1,10)

  res_abar <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar0, Anodes = Anodes)
  res_abar <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar1, Anodes = Anodes)
  checkException(
  res_abar <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar2, Anodes = Anodes)
  )
  res_abar <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar3, Anodes = Anodes)


  Anodes <- "A"%+%c(1:5)
  test.dat <- matrix(rnorm(50), nrow=10, ncol =5)
  abar0 <- 0
  abar1 <- c(0,1,1,0,1)
  abar2 <- matrix(c(0,1,1,0,1), nrow=nrow(test.dat), ncol = length(Anodes), byrow=TRUE)
  abar3 <- function(data, ...) 1
  abar4 <- function(data, ...) rep(0,5)
  abar5 <- function(data, ...) matrix(c(1,1,1,0,1), nrow=nrow(data), ncol = 5, byrow=FALSE)
  abar6 <- function(data, ...) {
    mat <- matrix(c(1,1,1,0,1), nrow=nrow(data), ncol = 5, byrow=FALSE)
    colnames(mat) <- "A"%+%c(1:5)
    return(mat)
  }
  abar7.a <- lapply(Anodes, function(Anode) return(def = function(data,...){c(0,1,1,0,1,0,1,1,0,1)}))
  abar7.b <- abar7.a; names(abar7.b) <- Anodes


  checkException(
    abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar0, Anodes = Anodes)
    )

  abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar1, Anodes = Anodes)
  checkEquals(is.matrix(abar_mat), TRUE)

  checkException(
  abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar2, Anodes = Anodes)
  )
  colnames(abar2) <- Anodes
  abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar2, Anodes = Anodes)
  colnames(abar2)[1] <- "A30"
  checkException(
  abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar2, Anodes = Anodes)
  )

  checkException(abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar3, Anodes = Anodes))

  checkException(abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar4, Anodes = Anodes))

  abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar6, Anodes = Anodes)

  checkException(abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar7.a, Anodes = Anodes))

  abar_mat <- tmlenet:::f.gen.A.star(data = test.dat, f.g_fun = abar7.b, Anodes = Anodes)
  checkEquals(is.matrix(abar_mat), TRUE)
}




