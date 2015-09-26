#***************************************************************************************
# TO DO: 
# (1) ADD AN EXAMPLE WITH STOCHASTIC INTERVENTION
# (2) ADD AN EXAMPLE WITH CONTINUOUS EXPOSURE

#***************************************************************************************
data(df_netKmax6) # Load the network data
Kmax <- 6 # Max number of friends in the network
#***************************************************************************************

#***************************************************************************************
# EXAMPLES OF INTERVENTION FUNCTIONS:
#***************************************************************************************
# Returns a function that will sample A with probability x:=P(A=1))
make_f.gstar <- function(x, ...) {
  eval(x)
  f.A_x <- function(data, ...){
    rbinom(n = nrow(data), size = 1, prob = x[1])
  }
  return(f.A_x)
}
# Deterministic f_gstar setting every A=0:
f.A_0 <- make_f.gstar(x = 0)
# Deterministic f_gstar setting every A=1:
f.A_1 <- make_f.gstar(x = 1)
# Stochastically sets (100*x)% of population to A=1 with probability 0.2:
f.A_.2 <- make_f.gstar(x = 0.2)

#***************************************************************************************
# DEFINING SUMMARY MEASURES:
#***************************************************************************************
def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
          def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
          def.sA(netA = A[[0:Kmax]])

#***************************************************************************************
# EVALUATING SUMMARY MEASURES FROM INPUT DATA
#***************************************************************************************
# A helper function that can pre-evaluate the summary measures on (O)bserved data:
res <- eval.summaries(sW = def_sW, sA = def_sA,  Kmax = 6, data = df_netKmax6,
  NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

#***************************************************************************************
# Specifying the clever covariate regressions hform.g0 and hform.gstar:
# Left side can consist of any summary names defined by def.sA (as linear terms)
# Right side can consist of any summary names defined by def.sW (as linear terms) & 'nF'
#***************************************************************************************
hform.g01 <- "netA ~ netW2 + sum.netW3 + nF"
hform.gstar1 <- "netA ~ sum.netW3"

# alternatives:
hform.g02 <- "netA + sum.netAW2 ~ netW2 + sum.netW3 + nF"
hform.g03 <- "sum.netAW2 ~ netW2 + sum.netW3"

#***************************************************************************************
# SPECIFYING the outcome regression Qform:
# Left side is ignored (with a warning if not equal to Ynode)
# Right side can by any summary names defined by def.sW, def.sA & 'nF'
#***************************************************************************************
Qform1 <- "Y ~ sum.netW3 + sum.netAW2"
Qform2 <- "Y ~ netA + netW + nF"
Qform3 <- "blah ~ netA + netW + nF"

#***************************************************************************************
# Estimate mean population outcome under deterministic intervention A=0 with 6 friends:
# ESTIMATION WITH regression formulas
#***************************************************************************************
# Note that Ynode is optional when Qform is specified;
options(tmlenet.verbose = FALSE) # set to TRUE to print status messages
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    Anode = "A", f_gstar1 = 0L,
                    IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

res_K6_1$EY_gstar1$estimates
res_K6_1$EY_gstar1$vars
res_K6_1$EY_gstar1$CIs

#***************************************************************************************
# Same as above but for covariate-based TMLE.
#***************************************************************************************
res_K6_2 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    Anode = "A", f_gstar1 = 0L,
                    IDnode = "IDs", NETIDnode = "Net_str", sep = ' ',
                    optPars = list(runTMLE = "tmle.covariate", n_MCsims = 10))

res_K6_2$EY_gstar1$estimates
res_K6_2$EY_gstar1$vars
res_K6_2$EY_gstar1$CIs

#***************************************************************************************
# SPECIFYING THE NETWORK AS A MATRIX OF FRIEND ROW NUMBERS:
#***************************************************************************************
net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = Kmax)

# generating the network matrix from input data:
NetInd_mat <- net_ind_obj$makeNetInd.fromIDs(Net_str = df_netKmax6[, "Net_str"],
                                              IDs_str = df_netKmax6[, "IDs"],
                                              sep = ' ')$NetInd
# number of friends:
nF <- net_ind_obj$nF

data(NetInd_mat_Kmax6)
all.equal(NetInd_mat, NetInd_mat_Kmax6) # TRUE

print(head(NetInd_mat))
print(head(nF))
print(all.equal(df_netKmax6[,"nFriends"], nF))

res_K6_net1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                        Qform = "Y ~ sum.netW3 + sum.netAW2",
                        hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                        hform.gstar = "netA ~ sum.netW3",
                        Anode = "A", f_gstar1 = f.A_0,
                        NETIDmat = NetInd_mat,
                        optPars = list(runTMLE = "tmle.intercept", n_MCsims = 10))

all.equal(res_K6_net1$EY_gstar1$estimates, res_K6_1$EY_gstar1$estimates)
all.equal(res_K6_net1$EY_gstar1$vars, res_K6_1$EY_gstar1$vars)
all.equal(res_K6_net1$EY_gstar1$CIs, res_K6_1$EY_gstar1$CIs)

#***************************************************************************************
# (*) EQUIVALENT WAYS OF SPECIFYING INTERVENTIONS f_gstar1/f_gstar2.
# (*) LOWERING THE DIMENSIONALITY OF THE SUMMARY MEASURES.
#***************************************************************************************
def_sW <- def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)
def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE)

# can define intervention by function f.A_0 that sets everyone's A to constant 0:
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = f.A_0,
                    IDnode = "IDs", NETIDnode = "Net_str")
res_K6_1$EY_gstar1$estimates

# equivalent way to define intervention f.A_0 is to just set f_gstar1 to 0:
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = 0L,
                    IDnode = "IDs", NETIDnode = "Net_str")
res_K6_1$EY_gstar1$estimates

# or set f_gstar1 to a vector of 0's of length nrow(data):
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = rep_len(0L, nrow(df_netKmax6)),
                    IDnode = "IDs", NETIDnode = "Net_str")
res_K6_1$EY_gstar1$estimates


#***************************************************************************************
# EXAMPLE WITH SIMULATED DATA FOR 2 FRIENDS AND 1 COVARIATE W1 (SIMULATION 1)
#***************************************************************************************
# data(df_netKmax2)
# head(df_netKmax2)
# #--------------------------------------------------------
# # Define regression formulas for Q and g
# #--------------------------------------------------------
# Qform <- "Y ~  W1 + A + netW1_1 + netW1_2 + netA_1 + netA_2 + nF"
# # gform <- "A ~  W1 + netW1_1 + netW1_2 + nF"

# #***************************************************************************************
# # Mean population outcome under stochastic intervention P(A=1)=0.2
# #***************************************************************************************
# Net_str <- df_netKmax6[, "Net_str"]
# IDs_str <- df_netKmax6[, "IDs"]
# net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = Kmax)
# net_ind_obj$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = ' ')
# NetInd_mat <- net_ind_obj$NetInd

# tmlenet_out2 <- tmlenet(data = df_netKmax2, Anode = "A", Ynode = "Y",
#                         Kmax = 2, IDnode = "IDs", NETIDnode = "Net_str",
#                         Qform = Qform, gform = gform,
#                         f.g1.star = f.A_.2,
#                         n_MCsims = 4000, n_samp_g0gstar = 100)

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

# #***************************************************************************************
# # Average treatment effect (ATE) for two interventions, f.g1.star: A=1 vs f.g2.star: A=0
# #***************************************************************************************
# tmlenet_out3 <- tmlenet(data = df_netKmax2, Anode = "A", Ynode = "Y",
#                         Kmax = 2, IDnode = "IDs", NETIDnode = "Net_str",
#                         Qform = Qform, gform = gform,
#                         f.g1.star = f.A_1, f.g2.star = f.A_0,
#                         n_MCsims = 4000, n_samp_g0gstar = 100)

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









