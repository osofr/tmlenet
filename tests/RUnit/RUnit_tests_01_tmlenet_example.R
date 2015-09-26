# ---------------------------------------------------------------------------------
# TEST SET 1
# SIMPLE NETWORK TMLE
# ---------------------------------------------------------------------------------

test.examples <- function() {
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

  #***************************************************************************************
  # # correct version(s):
  #***************************************************************************************
  # No Ynode:
 res_K6_1a <- tmlenet(data = df_netKmax6,
              Qform = "Y ~ sum.netW3 + sum.netAW2",
              hform = "netA ~ netW2 + sum.netW3 + nF",
              hform.gstar = "netA ~ sum.netW3",

              Anode = "A",

              Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
              f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

  # With Ynode:
  res_K6_1 <- tmlenet(data = df_netKmax6,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",

                    Anode = "A", Ynode = "Y",

                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept",n_MCsims = 10))

  tmle_idx <- rownames(res_K6_1$EY_gstar1$estimates)%in%"tmle"
  h_iptw_idx <- rownames(res_K6_1$EY_gstar1$estimates)%in%"h_iptw"
  gcomp_idx <- rownames(res_K6_1$EY_gstar1$estimates)%in%"gcomp"

  # Test estimates:
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[tmle_idx] - 0.5051903) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[h_iptw_idx] - 0.5065960) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[gcomp_idx] - 0.4970377) < 10^(-06))
  # Test asymptotic vars:
  checkTrue(abs(res_K6_1$EY_gstar1$vars[tmle_idx] - 0.0009268804) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$vars[h_iptw_idx] - 0.0021023317) < 10^(-06))
  # Test CIs:
  checkTrue((abs(res_K6_1$EY_gstar1$CIs[tmle_idx][1] - 0.4455197) < 10^(-06)) &
              (abs(res_K6_1$EY_gstar1$CIs[tmle_idx][2] - 0.5648608) < 10^(-06)))
  checkTrue((abs(res_K6_1$EY_gstar1$CIs[h_iptw_idx][1] - 0.4167293) < 10^(-06)) &
              (abs(res_K6_1$EY_gstar1$CIs[h_iptw_idx][2] - 0.5964627) < 10^(-06)))
  # res_K6_1$EY_gstar1$other.vars

  #----------------------------------------------------------------------------------
  # Example 2. Same as above but for covariate-based TMLE
  #----------------------------------------------------------------------------------
  res_K6_2 <- tmlenet(data = df_netKmax6,
                      Qform = "Y ~ sum.netW3 + sum.netAW2",
                      hform = "netA ~ netW2 + sum.netW3 + nF",
                      hform.gstar = "netA ~ sum.netW3",

                      Anode = "A", Ynode = "Y",
                      Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                      f_gstar1 = f.A_0, sW = def_sW, sA = def_sA,
                      optPars = list(runTMLE = "tmle.covariate",n_MCsims = 10))

  # Test estimates:
  checkTrue(abs(res_K6_2$EY_gstar1$estimates[tmle_idx] - 0.5053725) < 10^(-06))
  checkTrue(abs(res_K6_2$EY_gstar1$estimates[h_iptw_idx] - 0.5065960) < 10^(-06))
  checkTrue(abs(res_K6_2$EY_gstar1$estimates[gcomp_idx] - 0.4970377) < 10^(-06))
  # Test asymptotic vars:
  checkTrue(abs(res_K6_2$EY_gstar1$vars[tmle_idx] - 0.0009265246) < 10^(-06))
  checkTrue(abs(res_K6_2$EY_gstar1$vars[h_iptw_idx] - 0.0021023317) < 10^(-06))
  # Test CIs:
  checkTrue((abs(res_K6_2$EY_gstar1$CIs[tmle_idx][1] - 0.4457134) < 10^(-06)) &
              (abs(res_K6_2$EY_gstar1$CIs[tmle_idx][2] - 0.5650315) < 10^(-06)))
  checkTrue((abs(res_K6_2$EY_gstar1$CIs[h_iptw_idx][1] - 0.4167293) < 10^(-06)) &
              (abs(res_K6_2$EY_gstar1$CIs[h_iptw_idx][2] - 0.5964627) < 10^(-06)))
  # res_K6$EY_gstar1$other.vars

  # ================================================================
  # COMPARING OLD vs NEW OUTPUT
  # N=1,000
  # NEW INTERFACE FINAL RESULTS WITH MC EVAL (FULL MATCH TO OLD)
  # gIPTW and TMLE_gIPTW AREN'T IMPLEMENTED
  # ================================================================
  # with h method is subsetting (excluding degenerate outcome sA):
                               # old:   # new:
  # epsilon (covariate)      0.02549743 0.02549743
  # alpha (intercept)        0.05410938 0.05410938
  # iptw epsilon (covariate) 0.03556655 0.03556655

                # old:  # new:
  # tmle_A     0.5053725 0.5053725
  # tmle_B     0.5051903 0.5051903
  # iid.tmle_B 0.4475714 0.4475714
  # tmle_iptw  0.5123310 0.5123310
  # iptw_h     0.5065960 0.5065960
  # iptw       0.4910014 0.4910014
  # iid.iptw   0.4429414 0.4429414
  # mle        0.4970377 0.4970377

  #   fWi_init_A   fWi_init_B   fWi_star_A   fWi_star_B
  # -0.008334713 -0.008152518 -0.505372462 -0.505190268

  # NEW VARs:
  # gIPTW and TMLE_gIPTW AREN'T YET IMPLEMENTED
  # $EY_gstar1$vars
  #                      var
  # tmle_A      0.0009265246
  # tmle_B      0.0009268804
  # tmle_g_iptw 0.0069826701
  # h_iptw      0.0021023317
  # g_iptw      0.0000000000
  # mle         0.0000000000

  # NEW CIs
  # $EY_gstar1$CIs
  #             LBCI_0.025 UBCI_0.975
  # tmle_A       0.4457134  0.5650315
  # tmle_B       0.4455197  0.5648608
  # tmle_g_iptw -0.1637792  0.1637792
  # h_iptw       0.4167293  0.5964627
  # g_iptw       0.0000000  0.0000000
  # mle          0.4970377  0.4970377
  # > tmlenet_K6out2$EY_gstar1$other.vars
  #    var_iid.tmle_B var_tmleiptw_2ndO     var_iptw_2ndO var_tmle_A_Q.init var_tmle_B_Q.init
  #      0.0004350965      0.0846013137      0.0000000000      0.0008532711      0.0008535512

  #----------------------------------------------------------------------------------
  # Example 3. Same as Example 1 but with with true f_g0
  # *** Note that since f_g0 depends on (W1, netW1, netW2, netW3), these covariates also need to be added to sW summary measure ***
  #----------------------------------------------------------------------------------
  def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
              def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE) +
              def.sW(netW1 = W1[[0:Kmax]], netW3 = W3[[1:Kmax]])

  def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]]) * W2[[1:Kmax]]), replaceNAw0 = TRUE) +
              def.sA(netA = A[[0:Kmax]])

  res_K6_3 <- tmlenet(data = df_netKmax6,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept",f_g0 = f.A,n_MCsims = 10))

  # Test estimates:
  checkTrue(abs(res_K6_3$EY_gstar1$estimates[tmle_idx] - 0.5054745) < 10^(-06))
  checkTrue(abs(res_K6_3$EY_gstar1$estimates[h_iptw_idx] - 0.4668999) < 10^(-06))
  checkTrue(abs(res_K6_3$EY_gstar1$estimates[gcomp_idx] - 0.4970377) < 10^(-06))
  # Test asymptotic vars:
  checkTrue(abs(res_K6_3$EY_gstar1$vars[tmle_idx] - 0.0008937359) < 10^(-06))
  checkTrue(abs(res_K6_3$EY_gstar1$vars[h_iptw_idx] - 0.001738169) < 10^(-06))
  # Test CIs:
  checkTrue((abs(res_K6_3$EY_gstar1$CIs[tmle_idx][1] - 0.4468806) < 10^(-06)) &
            (abs(res_K6_3$EY_gstar1$CIs[tmle_idx][2] - 0.5640685) < 10^(-06)))
  checkTrue((abs(res_K6_3$EY_gstar1$CIs[h_iptw_idx][1] - 0.3851863) < 10^(-06)) &
            (abs(res_K6_3$EY_gstar1$CIs[h_iptw_idx][2] - 0.5486134) < 10^(-06)))
  # res_K6$EY_gstar1$other.vars

  #----------------------------------------------------------------------------------
  # Same as Example 1, but specifying the network with NETIDmat: a matrix of friend row numbers from the input data
  #----------------------------------------------------------------------------------
  Net_str <- df_netKmax6[, "Net_str"]
  IDs_str <- df_netKmax6[, "IDs"]
  net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = Kmax)
  net_ind_obj$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = ' ')
  NetInd_mat <- net_ind_obj$NetInd
  # save(NetInd_mat_Kmax6, file = "NetInd_mat_Kmax6.rda")

  nF <- net_ind_obj$nF
  print(head(NetInd_mat))
  print(head(nF))
  checkTrue(all.equal(df_netKmax6[,"nFriends"], nF))

  res_K6net <- tmlenet(data = df_netKmax6,
                      Qform = "Y ~ sum.netW3 + sum.netAW2",
                      hform = "netA ~ netW2 + sum.netW3 + nF",
                      hform.gstar = "netA ~ sum.netW3",

                      Anode = "A", Ynode = "Y",
                      Kmax = Kmax, NETIDmat = NetInd_mat,
                      f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept",n_MCsims = 10))

  checkTrue(all.equal(res_K6net$EY_gstar1$estimates, res_K6_1$EY_gstar1$estimates))
  checkTrue(all.equal(res_K6net$EY_gstar1$vars, res_K6_1$EY_gstar1$vars))
  checkTrue(all.equal(res_K6net$EY_gstar1$CIs, res_K6_1$EY_gstar1$CIs))
  checkTrue(all.equal(res_K6net$EY_gstar1$other.vars, res_K6_1$EY_gstar1$other.vars))
}

#***************************************************************************************
# EXAMPLE WITH SIMULATED DATA FOR 2 FRIENDS AND 1 COVARIATE W1 (SIMULATION 1)
#***************************************************************************************
# data(sample_network_k2)
# load(file="./sample_network_k2.RData")
# head(sample_network_k2)

#--------------------------------------------------------
# Define regression formulas for Q and g
# ****IMPORTANT****: 
#	use notation netVAR_1 to refer to covariate VAR of the 1st friend
# 	netVAR_2 to refer to covariate VAR of the 2nd friend and so on...
#--------------------------------------------------------
# Qform <- "Y ~  W1 + A + netW1_1 + netW1_2 + netA_1 + netA_2 + nF"
# gform <- "A ~  W1 + netW1_1 + netW1_2 + nF"

#----------------------------------------------------------------------------------
# Example 2. Mean population outcome under stochastic intervention P(A=1)=0.2
# OLD. REMOVE OR MODIFY.
#----------------------------------------------------------------------------------
# tmlenet_out2 <- tmlenet(data=sample_network_k2, Anode="A", Ynode="Y",
# 						Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
# 						f.g1.star=f.A_x, f.g1_args=list(x=0.2),
# 						n_MCsims=4000, n_samp_g0gstar=100)

# # TMLE estimate and iid IC-based 95% CI:
# tmlenet_out2$estimates$EY_g1.star$tmle_B
# tmlenet_out2$estimates$EY_g1.star$CI_tmle_B_iidIC

# # Efficient IPTW (h) and iid IC-based 95% CI:
# tmlenet_out2$estimates$EY_g1.star$iptw_h
# tmlenet_out2$estimates$EY_g1.star$CI_iptw_h_iidIC

# # Inefficient IPTW (g) + (two CIs, less conservative and more conservative)
# tmlenet_out2$estimates$EY_g1.star$iptw
# tmlenet_out2$estimates$EY_g1.star$CI_iptw_iidIC_1stO
# tmlenet_out2$estimates$EY_g1.star$CI_iptw_iidIC_2ndO

# # MLE
# tmlenet_out2$estimates$EY_g1.star$mle

#----------------------------------------------------------------------------------
# Example 3. Average treatment effect (ATE) for two interventions, f.g1.star: A=1 vs f.g2.star: A=0
# OLD. REMOVE OR MODIFY.
#----------------------------------------------------------------------------------
# tmlenet_out3 <- tmlenet(data=sample_network_k2, Anode="A", Ynode="Y",
# 						Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
# 						f.g1.star=f.A_1, f.g1_args=NULL, f.g2.star=f.A_0, f.g2_args=NULL,
# 						n_MCsims=4000, n_samp_g0gstar=100)

# # TMLE estimate for ATE + 95% CI
# tmlenet_out3$estimates$ATE$tmle_B
# tmlenet_out3$estimates$ATE$CI_tmle_B_iidIC

# # Efficient IPTW (h) and iid IC-based 95% CI:
# tmlenet_out3$estimates$ATE$iptw_h
# tmlenet_out3$estimates$ATE$CI_iptw_h_iidIC

# # Inefficient IPTW (g) + (two CIs, less conservative and more conservative)
# tmlenet_out3$estimates$ATE$iptw
# tmlenet_out3$estimates$ATE$CI_iptw_iidIC_1stO
# tmlenet_out3$estimates$ATE$CI_iptw_iidIC_2ndO

# # MLE
# tmlenet_out3$estimates$ATE$mle
