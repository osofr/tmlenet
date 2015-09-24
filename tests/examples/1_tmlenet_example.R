#***************************************************************************************
data(df_netKmax6) # Load the network data
Kmax <- 6 # Max number of friends in the network
#***************************************************************************************

#***************************************************************************************
# Example 1. 
# Mean population outcome under deterministic intervention A=0 with 6 friends.
# Intercept based TMLE.
#***************************************************************************************

# *************
# TO BE REMOVED (no longer need Wnodes argument in tmlenet())
# *************
Wnodes <- c("W1", "W2", "W3")

# SUMMARY MEASURES:
def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
  def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
          def.sA(netA = A[[0:Kmax]])

# POSSIBLE INTERVENTION FUNCTIONS:
# Set x% of community to A=1 (returns A sampled with probability P(A=1)):
# f.A_x <- function(data, x, ...) rep(x, nrow(data))
f.A_x <- function(data, x, ...) rbinom(n = nrow(data), size = 1, prob = x[1])
# Deterministically set every A=0:
f.A_0 <- function(data, ...) f.A_x(data, 0, ...)
# Deterministically set every A=1:
f.A_1 <- function(data, ...) f.A_x(data, 1, ...)

# alternative ways to pass summary measures:
# sW = list("W1[[0]]", "W2[[0:Kmax]]", "W3[[0:Kmax]]", 
            # netW1_sum = "rowSums(W1[[1:Kmax]]"), 
            # netW2_sum = "rowSums(W2[[1:Kmax]])", 
            # netW3_sum = "rowSums(W3[[1:Kmax]])"), 
# sA = list("A[[0:Kmax]]", 
            # sum_1mAW2_nets = "rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))")

# NOT IMPLEMENTED YET. 
# A helper function that can pre-evaluate the summary measures on (O)bserved data 
# (data.frame)
# This will help when examining the data and playing with various summary measures, 
# prior to running the tmletnet() function
# res <- eval.summaries(summaries = def_sA, Odata = df_netKmax6, Kmax = Kmax, 
# NETIDnode = "Net_str", IDnode = "IDs")

# --------------------------------------------------
# # NEW INTERFACE FOR SPECIFYING hform, Qform, gform allows including the summary 
# measure names
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

options(tmlenet.verbose = FALSE)
res_K6_1 <- tmlenet(data = df_netKmax6, Anode = "A", Wnodes = Wnodes, Ynode = "Y",
                  Kmax = Kmax,
                  IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                  f_gstar1 = f.A_0,
                  sW = def_sW, sA = def_sA,
                  Qform = "Y ~ sum.netW3 + sum.netAW2",
                  hform = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",
                  optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))


res_K6_1$EY_gstar1$estimates
res_K6_1$EY_gstar1$vars
res_K6_1$EY_gstar1$CIs
res_K6_1$EY_gstar1$other.vars




#***************************************************************************************
# Example 2. 
# Same as above but for covariate-based TMLE.
#***************************************************************************************
res_K6_2 <- tmlenet(data = df_netKmax6, Anode = "A", Wnodes = Wnodes, Ynode = "Y",
                  Kmax = Kmax,
                  IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                  f_gstar1 = f.A_0,
                  sW = def_sW, sA = def_sA,
                  Qform = "Y ~ netW3_sum + sum_1mAW2_nets",
                  hform = "netA ~ netW2 + netW3_sum + nF",
                  hform.gstar = "netA ~ netW3_sum",
                  optPars = list(runTMLE = "tmle.covariate", n_MCsims = 10))
res_K6_2$EY_gstar1$estimates
res_K6_2$EY_gstar1$vars
res_K6_2$EY_gstar1$CIs
res_K6_2$EY_gstar1$other.vars

#***************************************************************************************
# Example 3. 
# Same as Example 1, but specifying the network with a matrix of friend row numbers.
#***************************************************************************************
Net_str <- df_netKmax6[, "Net_str"]
IDs_str <- df_netKmax6[, "IDs"]
net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = Kmax)
net_ind_obj$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = ' ')
NetInd_mat <- net_ind_obj$NetInd

data(NetInd_mat_Kmax6)
all.equal(NetInd_mat, NetInd_mat_Kmax6) # TRUE

nF <- net_ind_obj$nF
print(head(NetInd_mat))
print(head(nF))
print(all.equal(df_netKmax6[,"nFriends"], nF))

res_K6_net1 <- tmlenet(data = df_netKmax6, Anode = "A", Wnodes = Wnodes,
                    Ynode = "Y",
                    Kmax = Kmax,
                    NETIDmat = NetInd_mat,
                    f_gstar1 = f.A_0,
                    sW = def_sW, sA = def_sA,
                    Qform = "Y ~ netW3_sum + sum_1mAW2_nets",
                    hform = "netA ~ netW2 + netW3_sum + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

all.equal(res_K6_net1$EY_gstar1$estimates, res_K6_1$EY_gstar1$estimates)
all.equal(res_K6_net1$EY_gstar1$vars, res_K6_1$EY_gstar1$vars)
all.equal(res_K6_net1$EY_gstar1$CIs, res_K6_1$EY_gstar1$CIs)
all.equal(res_K6_net1$EY_gstar1$other.vars, res_K6_1$EY_gstar1$other.vars)

#***************************************************************************************
# Example 2. Mean population outcome under deterministic intervention A=1 with 6 friends
# OLD. REMOVE OR MODIFY.
#***************************************************************************************
# tmlenet_K6out2 <- tmlenet(data=df_netKmax6, Anode="A", Wnodes=Wnodes, Ynode="Y", 
# Kmax=Kmax, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform, h_form=hform,
# f.g1.star=f.A_1, f.g1_args=NULL, n_MCsims=10, n_samp_g0gstar=10)

# tmlenet_K6out2$estimates$EY_g1.star$tmle_B
# tmlenet_K6out2$estimates$EY_g1.star$CI_tmle_B_iidIC


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

#***************************************************************************************
# Example 1. Mean population outcome under deterministic intervention A=0
# OLD. REMOVE OR MODIFY.
#***************************************************************************************
# tmlenet_out1 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", 
# Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform, 
# f.g1.star=f.A_0, f.g1_args=NULL)

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

#***************************************************************************************
# Example 2. Mean population outcome under stochastic intervention P(A=1)=0.2
# OLD. REMOVE OR MODIFY.
#***************************************************************************************
# tmlenet_out2 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", 
# Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
# f.g1.star=f.A_x, f.g1_args=list(x=0.2),
# n_MCsims=4000, n_samp_g0gstar=100)

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

#***************************************************************************************
# Example 3. Average treatment effect (ATE) for two interventions, 
# f.g1.star: A=1 vs f.g2.star: A=0
# OLD. REMOVE OR MODIFY.
#***************************************************************************************
# tmlenet_out3 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", 
# Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
# f.g1.star=f.A_1, f.g1_args=NULL, f.g2.star=f.A_0, f.g2_args=NULL,
# n_MCsims=4000, n_samp_g0gstar=100)

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
