#***************************************************************************************
# TO DO: 
# (1) ADD AN EXAMPLE WITH STOCHASTIC INTERVENTION
# (2) ADD AN EXAMPLE WITH CONTINUOUS EXPOSURE
#***************************************************************************************
data(df_netKmax6) # Load the network data
Kmax <- 6 # Max number of friends in the network

#***************************************************************************************
# Example 1. 
# Mean population outcome under deterministic intervention A=0 with 6 friends.
# Intercept based TMLE.
#***************************************************************************************

#***************************************************************************************
# POSSIBLE INTERVENTION FUNCTIONS:
#***************************************************************************************
# Set x% of community to A=1 (returns A sampled with probability P(A=1)):
f.A_x <- function(data, x, ...) rbinom(n = nrow(data), size = 1, prob = x[1])
# Deterministically set every A=0:
f.A_0 <- function(data, ...) f.A_x(data, 0, ...)
# Deterministically set every A=1:
f.A_1 <- function(data, ...) f.A_x(data, 1, ...)

#***************************************************************************************
# SUMMARY MEASURES:
#***************************************************************************************
def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
  def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
          def.sA(netA = A[[0:Kmax]])

#***************************************************************************************
# EVALUATING SUMMARY MEASURES FOR INPUT DATA
# A helper function that can pre-evaluate the summary measures on (O)bserved data:
#***************************************************************************************
res <- eval.summaries(sW = def_sW, sA = def_sA,  Kmax = 6, data = df_netKmax6,
  NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

#***************************************************************************************
# ESTIMATION
#***************************************************************************************
options(tmlenet.verbose = FALSE)
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                  Anode = "A", Ynode = "Y", f_gstar1 = f.A_0,
                  IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                  optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

res_K6_1$EY_gstar1$estimates
res_K6_1$EY_gstar1$vars
res_K6_1$EY_gstar1$CIs
res_K6_1$EY_gstar1$other.vars

#***************************************************************************************
# SPECIFYING the clever covariate regressions hform and hform.gstar:
# Left side can consist of any summary names defined by def.sA (as linear terms)
# Right side can consist of any summary names defined by def.sW (as linear terms) & 'nF'
#***************************************************************************************
hform1 <- "netA ~ netW2 + sum.netW3 + nF"
hform.gstar1 <- "netA ~ sum.netW3"

# alternatives:
hform2 <- "netA + sum.netAW2 ~ netW2 + sum.netW3 + nF"
hform3 <- "sum.netAW2 ~ netW2 + sum.netW3"

#***************************************************************************************
# SPECIFYING the outcome regression Qform:
# Left side is ignored (with a warning if not equal to Ynode)
# Right side can by any summary names defined by def.sW, def.sA & 'nF'
#***************************************************************************************
Qform1 <- "Y ~ sum.netW3 + sum.netAW2"
Qform2 <- "Y ~ netA + netW + nF"
Qform3 <- "blah ~ netA + netW + nF"

#***************************************************************************************
# ESTIMATION WITH regression formulas:
# Note that Ynode is optional when Qform is specified;
#***************************************************************************************
options(tmlenet.verbose = FALSE)
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                  Qform = "Y ~ sum.netW3 + sum.netAW2",
                  hform = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",
                  Anode = "A",
                  f_gstar1 = f.A_0,
                  IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                  optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

res_K6_1$EY_gstar1$estimates
res_K6_1$EY_gstar1$vars
res_K6_1$EY_gstar1$CIs
res_K6_1$EY_gstar1$other.vars

#***************************************************************************************
# Example 2. 
# Same as above but for covariate-based TMLE.
#***************************************************************************************
res_K6_2 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                  Qform = "Y ~ sum.netW3 + sum.netAW2",
                  hform = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",
                  Anode = "A", Ynode = "Y", f_gstar1 = f.A_0,
                  IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
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

res_K6_net1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    Anode = "A", Ynode = "Y", f_gstar1 = f.A_0,
                    NETIDmat = NetInd_mat,
                    optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

all.equal(res_K6_net1$EY_gstar1$estimates, res_K6_1$EY_gstar1$estimates)
all.equal(res_K6_net1$EY_gstar1$vars, res_K6_1$EY_gstar1$vars)
all.equal(res_K6_net1$EY_gstar1$CIs, res_K6_1$EY_gstar1$CIs)
all.equal(res_K6_net1$EY_gstar1$other.vars, res_K6_1$EY_gstar1$other.vars)

#***************************************************************************************
# EXAMPLE WITH SIMULATED DATA FOR 2 FRIENDS AND 1 COVARIATE W1 (SIMULATION 1)
#***************************************************************************************
# data(sample_network_k2)
# load(file="./sample_network_k2.RData")
# head(sample_network_k2)

#--------------------------------------------------------
# Define regression formulas for Q and g
# ****IMPORTANT****: 
# use notation netVAR_1 to refer to covariate VAR of the 1st friend
#   netVAR_2 to refer to covariate VAR of the 2nd friend and so on...
#--------------------------------------------------------
# Qform <- "Y ~  W1 + A + netW1_1 + netW1_2 + netA_1 + netA_2 + nF"
# gform <- "A ~  W1 + netW1_1 + netW1_2 + nF"

#***************************************************************************************
# Mean population outcome under stochastic intervention P(A=1)=0.2
#***************************************************************************************
# tmlenet_out2 <- tmlenet(data=sample_network_k2, Anode="A", Ynode="Y",
#             Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
#             f.g1.star=f.A_x, f.g1_args=list(x=0.2),
#             n_MCsims=4000, n_samp_g0gstar=100)

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
# Average treatment effect (ATE) for two interventions, f.g1.star: A=1 vs f.g2.star: A=0
#***************************************************************************************
# tmlenet_out3 <- tmlenet(data=sample_network_k2, Anode="A", Ynode="Y",
#             Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
#             f.g1.star=f.A_1, f.g1_args=NULL, f.g2.star=f.A_0, f.g2_args=NULL,
#             n_MCsims=4000, n_samp_g0gstar=100)

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
