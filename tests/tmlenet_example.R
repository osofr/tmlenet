# rm(list=ls())

#--------------------------------------------------------
# NOTE: argument n_MCsims specifies the number of Monte-Carlo sims tmlenet needs to run. If running time is too slow try lower that number.
#--------------------------------------------------------

# Set x% of community to A=1 (returns probability P(A=1))
f.A_x <- function(data, x, ...) rep(x, nrow(data))
# Deterministically set every A=0
f.A_0 <- function(data, ...) f.A_x(data, 0, ...)
# Deterministically set every A=1
f.A_1 <- function(data, ...) f.A_x(data, 1, ...)

#***************************************************************************************
# EXAMPLE WITH SIMULATED DATA FOR 6 FRIENDS AND 3 W's (SIMULATION 3)
#***************************************************************************************
# library(tmle)
# library(locfit)
# library(xtable)
library(bigmemory)
library(biganalytics)
library(plyr)
options(bigmemory.typecast.warning = FALSE)

kmax <- 6	# Max number of friends in the network?
# simulate a dataset first
source("../datgen_nets/sim3_datgen_k6.R")
set.seed(543)
n <- 1000
# df_K6 <-gendata_pop(nC=1, n_arr=1000, k_arr=kmax, EC_arr=EC, f.g_list="f.A", f.g_args_list=list(NULL))
t <- system.time(df_K6 <- gendata_pop(nC=1, n_arr=n, k_arr=kmax, EC_arr=EC, f.g_list="f.A", f.g_args_list=list(NULL)))
t
# for n=10K:
#   user  system elapsed 
# 33.117   0.976  33.898 

head(df_K6)
  # IDs Y nFriends W1 W2 W3 netW1_sum netW2_sum netW3_sum A                  Net_str
# 1  I1 1        1  3  0  1         2         1         0 1                     I537
# 2  I2 1        5  3  1  0        16         4         3 0    I6 I58 I595 I641 I654
# 3  I3 0        5  5  0  0        18         3         3 0 I163 I637 I650 I722 I783
# 4  I4 1        2  2  1  0         6         1         0 1                 I49 I995
# 5  I5 1        3  3  1  1        11         3         1 0           I358 I369 I762
# 6  I6 1        2  2  0  1         6         2         0 1                I682 I917

class(df_K6$A) # [1] "integer"
class(df_K6$nFriends) # [1] "numeric"
table(df_K6$W1)
 #   0    1    2    3    4    5 
 # 475 1772 2900 2615 1718  520 
c(mean(df_K6$W1), mean(df_K6$W2), mean(df_K6$W3))
# [1] 2.4889 0.5719 0.6031

mean(df_K6$A) # [1] 0.198
mean(df_K6$Y) # [1] 0.435
mean(df_K6$nFriends) # [1] 3.307

# --------------------------------------------------
# # NEW INTERFACE FOR SPECIFYING hform, Qform, gform allows including the summary measure names
# --------------------------------------------------
# testform1 <- as.formula("sA + sA2 ~ sW1 + netW3_sum")
# testform2 <- as.formula("netA ~ netW2 + netW3_sum")
# testform <- testform1
# testterms <- terms(testform)
# # Getting predictor sW names:
# (sW.names <- attributes(testterms)$term.labels)
# sW.names.alt <- colnames(attributes(testterms)$factors)
# assert_that(all(sW.names == sW.names.alt))
# # Getting outcome sA names:
# (out.var <- rownames(attributes(testterms)$factors)[1]) # character string
# out.vars.form <- as.formula(". ~ " %+% out.var)
# out.vars.terms <- terms(out.vars.form)
# (sA.names <- attributes(out.vars.terms)$term.labels)

# --------------------------------------------------
# SUMMARY MEASURES
# --------------------------------------------------
# NOT IMPLEMENTED YET. 
# A helper function that can pre-evaluate the summary measures on (O)bserved data (data.frame)
# This will help when examining the data and playing with various summary measures, prior to running the tmletnet() function
# res <- eval.summaries(summaries = def_sA, Odata = df_K6, Kmax = kmax, NETIDnode = "Net_str", IDnode = "IDs")
# --------------------------------------------------

#----------------------------------------------------------------------------------
# Example 1. Mean population outcome under deterministic intervention A=0 with 6 friends
#----------------------------------------------------------------------------------
Wnodes <- c("W1", "W2", "W3", "netW1_sum", "netW2_sum", "netW3_sum")
head(df_K6)

def_sW <- def.sW(netW2 = W2[[1:Kmax]], noname = TRUE) + 
            def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceNAw0 = TRUE)
            
def_sA <- def.sA(sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]), replaceNAw0 = TRUE) +
            def.sA(netA = A[[0:Kmax]], noname = TRUE)

(Qform.depr <- "Y ~ I(" %+% paste(str_c(netvar("W2", (1:6)), "*", "(1-",netvar("A", (1:6)), ")"), collapse = "+") %+% ")" %+% "+"  %+% "netW3_sum")
(gform.depr <- "A ~  W1 + netW1_sum + netW2_sum + netW3_sum + nFriends")
(hform.depr <- "sA ~ " %+% paste(netvar("W2", (1:6)), collapse = "+") %+% " + netW3_sum + nFriends")

system.time(
tmlenet_K6out2 <- tmlenet(data = df_K6, Anode = "A", Wnodes = Wnodes, Ynode = "Y", nFnode = "nFriends",
                          Kmax = kmax, 
                          IDnode = "IDs", 
                          NETIDnode = "Net_str",
                          # NETIDnode = NULL,
                          f_gstar1 = f.A_0,

                          # OLD regs (TO BE REMOVED):
                          Qform.depr = Qform.depr, hform.depr = hform.depr, #gform.depr = gform.depr,  # remove

                          sW = def_sW, sA = def_sA,
                          # new way to specify regressions:
                          Qform = "Y ~ netW3_sum + sum_1mAW2_nets",
                          hform = "netA ~ netW2 + netW3_sum + nFriends",
                          hform.gstar = "netA ~ netW3_sum",
                          gform = "A ~  W1 + netW1_sum + netW2_sum + netW3_sum + nFriends",
                          opt.params = list(
                            onlyTMLE_B = FALSE,  # remove
                            # f_g0 = f.A, # tested, works
                            n_MCsims = 10
                          )
                          ))
                          # alternative way to pass summary measures:
                          # sW = list("W1[[0]]", "W2[[0:Kmax]]", "W3[[0:Kmax]]", netW1_sum = "rowSums(W1[[1:Kmax]]"), netW2_sum = "rowSums(W2[[1:Kmax]])", netW3_sum = "rowSums(W3[[1:Kmax]])"), 
                          # sA = list("A[[0:Kmax]]", sum_1mAW2_nets = "rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))")

tmlenet_K6out2$EY_gstar1$estimates
tmlenet_K6out2$EY_gstar1$vars
tmlenet_K6out2$EY_gstar1$CIs
tmlenet_K6out2$EY_gstar1$other.vars

# ================================================================
# COMPARING OLD vs NEW OUTPUT
# N=1,000
# ================================================================
# with h method is subsetting (excluding degenerate outcome sA):
                             # old: 	# new:
# epsilon (covariate)      0.02549743 0.02549743
# alpha (intercept)        0.05410938 0.05410938
# iptw epsilon (covariate) 0.03556655 0.03556655

              # old:	# new:
# tmle_A     0.5053725 0.5053725
# tmle_B     0.5051903 0.5051903
# iid.tmle_B 0.4475714 0.4475714
# tmle_iptw  0.5123310 0.5123310
# iptw_h     0.5065960 0.5065960
# iptw       0.4910014 0.4910014
# iid.iptw   0.4429414 0.4429414
# mle        0.4970377 0.4970377
# ================================================================
# NEW INTERFACE FINAL RESULTS WITH MC EVAL (FULL MATCH TO OLD)
# gIPTW and TMLE_gIPTW AREN'T YET IMPLEMENTED
# ================================================================
#                         h_iptw
# epsilon (covariate) 0.02549743
# alpha (intercept)   0.05410938
# [1] "time to run MCS: "
#    user  system elapsed 
#   0.099   0.010   0.109 
#   fWi_init_A   fWi_init_B   fWi_star_A   fWi_star_B 
# -0.008334713 -0.008152518 -0.505372462 -0.505190268 
# [1] "new MC.ests vec: "
#      tmle_A      tmle_B tmle_g_iptw      h_iptw      g_iptw         mle 
#   0.5053725   0.5051903   0.0000000   0.5065960   0.0000000   0.4970377 
# [1] "new MC.ests mat: "
#              estimate
# tmle_A      0.5053725
# tmle_B      0.5051903
# tmle_g_iptw 0.0000000
# h_iptw      0.5065960
# g_iptw      0.0000000
# mle         0.4970377

# ================================================================
# COMPARING OLD vs NEW Vars & CIs:
# ================================================================

# OLD tmle results object:
# > tmlenet_K6out2$EY_gstar1$vars
#                      var
# tmle_A      0.0009265246
# tmle_B      0.0009268804
# tmle_g_iptw 0.0005006418
# h_iptw      0.0021023317
# g_iptw      0.0060034386
# mle         0.0000000000
# > tmlenet_K6out2$EY_gstar1$CIs
#             LBCI_0.025 UBCI_0.975
# tmle_A       0.4457134  0.5650315
# tmle_B       0.4455197  0.5648608
# tmle_g_iptw  0.4698876  0.5575961
# h_iptw       0.4167293  0.5964627
# g_iptw       0.3669720  0.6706953
# mle          0.4970377  0.4970377
# > tmlenet_K6out2$EY_gstar1$other.vars
# var_iid.tmle_B var_tmleiptw_2ndO     var_iptw_2ndO var_tmle_A_Q.init var_tmle_B_Q.init 
#    0.0004350965      0.0003162110      0.0789868609      0.0008532711      0.0008535512 

# NEW CIs:
# gIPTW and TMLE_gIPTW AREN'T YET IMPLEMENTED
# $EY_gstar1$vars
#                      var
# tmle_A      0.0009265246
# tmle_B      0.0009268804
# tmle_g_iptw 0.0069826701
# h_iptw      0.0021023317
# g_iptw      0.0000000000
# mle         0.0000000000
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
# Example 2. Mean population outcome under deterministic intervention A=1 with 6 friends
# OLD. REMOVE OR MODIFY.
#----------------------------------------------------------------------------------
# tmlenet_K6out2 <- tmlenet(data=df_K6, Anode="A", Wnodes=Wnodes, Ynode="Y", nFnode="nFriends",
# 						Kmax=kmax, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform, h_form=hform,
# 						f.g1.star=f.A_1, f.g1_args=NULL, n_MCsims=10, n_samp_g0gstar=10)

# tmlenet_K6out2$estimates$EY_g1.star$tmle_B
# tmlenet_K6out2$estimates$EY_g1.star$CI_tmle_B_iidIC


#***************************************************************************************
# EXAMPLE WITH SIMULATED DATA FOR 2 FRIENDS AND 1 COVARIATE W1 (SIMULATION 1)
#***************************************************************************************
# data(sample_network_k2)
load(file="./sample_network_k2.RData")
head(sample_network_k2)

#--------------------------------------------------------
# Define regression formulas for Q and g
# ****IMPORTANT****: 
#	use notation netVAR_1 to refer to covariate VAR of the 1st friend
# 	netVAR_2 to refer to covariate VAR of the 2nd friend and so on...
#--------------------------------------------------------
Qform <- "Y ~  W1 + A + netW1_1 + netW1_2 + netA_1 + netA_2 + nFriends"
gform <- "A ~  W1 + netW1_1 + netW1_2 + nFriends"

#----------------------------------------------------------------------------------
# Example 1. Mean population outcome under deterministic intervention A=0
# OLD. REMOVE OR MODIFY.
#----------------------------------------------------------------------------------
# tmlenet_out1 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", nFnode="nFriends",
# 						Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform, 
# 						f.g1.star=f.A_0, f.g1_args=NULL)

# # TMLE estimate and iid IC-based 95% CI:
# tmlenet_out1$estimates$EY_g1.star$tmle_B
# tmlenet_out1$estimates$EY_g1.star$CI_tmle_B_iidIC

# # Efficient IPTW (h) and iid IC-based 95% CI:
# tmlenet_out1$estimates$EY_g1.star$iptw_h
# tmlenet_out1$estimates$EY_g1.star$CI_iptw_h_iidIC

# # Inefficient IPTW (g) + (two CIs, less conservative and more conservative)
# tmlenet_out1$estimates$EY_g1.star$iptw
# tmlenet_out1$estimates$EY_g1.star$CI_iptw_iidIC_1stO
# tmlenet_out1$estimates$EY_g1.star$CI_iptw_iidIC_2ndO

# # MLE
# tmlenet_out1$estimates$EY_g1.star$mle

#----------------------------------------------------------------------------------
# Example 2. Mean population outcome under stochastic intervention P(A=1)=0.2
# OLD. REMOVE OR MODIFY.
#----------------------------------------------------------------------------------
# tmlenet_out2 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", nFnode="nFriends",
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
# tmlenet_out3 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", nFnode="nFriends",
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