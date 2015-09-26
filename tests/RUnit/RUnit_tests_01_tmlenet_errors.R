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

  #***************************************************************************************
  # EXAMPLE WITH SIMULATED DATA FOR 6 FRIENDS AND 3 W's (OLD SIMULATION 3)
  #***************************************************************************************
  tmlenet:::checkpkgs(pkgs = c("stringr"))
  require(stringr)

  # Max number of friends in the network:
  Kmax <- 6
  # Load simulation function:
  source("./datgen_nets/sim3_datgen_k6.R")  # to load from inside run-time test dir
  # source("../datgen_nets/sim3_datgen_k6.R") # to load from current file dir
  # Simulate network data:
  # set.seed(543)
  n <- 1000
  df_netKmax6 <- gendata_pop(nC = 1, n_arr = n, k_arr = Kmax, EC_arr = EC, f.g_list = "f.A", f.g_args_list = list(NULL), rndseed = 543)
  # save(df_netKmax6, file = "./df_netKmax6.rda")

  print(head(df_netKmax6))
  #   IDs W1 W2 W3 A Y nFriends                  Net_str
  # 1  I1  3  0  1 1 1        1                     I537
  # 2  I2  3  1  0 0 1        5    I6 I58 I595 I641 I654
  # 3  I3  5  0  0 0 0        5 I163 I637 I650 I722 I783
  # 4  I4  2  1  0 1 1        2                 I49 I995
  # 5  I5  3  1  1 0 1        3           I358 I369 I762
  # 6  I6  2  0  1 1 1        2                I682 I917
  print(class(df_netKmax6$A)) # [1] "integer"
  print(class(df_netKmax6$nFriends)) # [1] "numeric"
  print(table(df_netKmax6$W1))
   #  0   1   2   3   4   5
   # 50 170 302 242 180  56
  print(c(mean(df_netKmax6$W1), mean(df_netKmax6$W2), mean(df_netKmax6$W3)))
  # [1] 2.500 0.553 0.589
  print(mean(df_netKmax6$A)) # [1] 0.198
  print(mean(df_netKmax6$Y)) # [1] 0.435
  print(mean(df_netKmax6$nFriends)) # [1] 3.307

  #----------------------------------------------------------------------------------
  # Example 1. Mean population outcome under deterministic intervention A=0 with 6 friends
  # Intercept based TMLE
  #----------------------------------------------------------------------------------
  options(tmlenet.verbose = FALSE)

  def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
            def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

  def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
            def.sA(netA = A[[0:Kmax]])

  # Test for non-existing predictors in hform/hform.gstar:
  checkException(
    res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform = "netA ~ netW2 + netW3_sum + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )
  # Test for non-existing predictors in hform/hform.gstar:
  checkException(
    res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )
  # Test for non-existing outcomes in hform/hform.gstar:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform = "sum.netW3 ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netW3 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )
  # Test for different outcomes in hform/hform.gstar (non-existing for hform.gstar):
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netW3 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )
  # Test for existing but different outcomes in hform/hform.gstar:
  checkException(
       res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netAW2 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )
  # Throws exception since netW3_sum, sum_1mAW2_nets from Qform don't exist:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Qform = " blah ~ netW3_sum + sum_1mAW2_nets",
                  hform = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",

                  Anode = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                  f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )
  # Throws exception since Ynode arg is omitted, but blah from Qform LHS doesn't exist:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Qform = " blah ~ sum.netW3 + sum.netAW2",
                  hform = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",

                  Anode = "A",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                  f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))
  )