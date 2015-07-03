rm(list=ls())

#------------------------------------
# Packages / package functions used so far inside tmlenet.
# If you are using just a few functions from another package,
# my recommendation is to note the package name in the Imports: field of the DESCRIPTION
# file and call the function(s) explicitly using ::, e.g., pkg::fun().

# If you are using functions repeatedly, you can avoid :: by importing the function with @importFrom pgk fun. 
# This also has a small performance benefit, because :: adds approximately 5 µs to function evaluation time.

# Alternatively, if you are repeatedly using many functions from another package, you can import all of them using @import package.
# This is the least recommended solution because it makes your code harder to read (you can’t tell where a function is coming from), 
# and if you @import many packages, it increases the chance of conflicting function names.

#-----------------------------------
# ADD tmlenet SOURCE CODE
library(R6)
library(assertthat)
# library(stringr)
# library(speedglm)
options(echo = TRUE)
options(width = 160)

library(devtools)
# load all R files in /R and datasets in /data. Ignores NAMESPACE:
load_all("../")
# source("BinarySummary_class.R")
# source("DatNet_class.R")
# source("Summaries_class.R")
# source("summarymeasures.R")
# source("tmlenet.R")

#--------------------------------------------------------
# NOTE: arguments n_MCsims and n_samp_g0gstar specify the number of Monte-Carlo sims tmlenet needs to run
# if things are taking too long try lower numbers (esspecially for n_MCsims)
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
library(locfit)
library(xtable)
library(bigmemory)
library(biganalytics)
library(plyr)
options(bigmemory.typecast.warning=FALSE)

kmax <- 6	# Max # of friends (K)?
# simulate a dataset first
source("./datgen_nets/sim3_datgen_k6.R")
set.seed(543)
# df_K6 <-gendata_pop(nC=1, n_arr=100, k_arr=kmax, EC_arr=EC, f.g_list="f.A", f.g_args_list=list(NULL))
df_K6 <- gendata_pop(nC=1, n_arr=1000, k_arr=kmax, EC_arr=EC, f.g_list="f.A", f.g_args_list=list(NULL))
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
mean(df_K6$A) # [1] 0.198
mean(df_K6$Y) # [1] 0.435
mean(df_K6$nFriends) # [1] 3.307

Wnodes <- c("W1", "W2", "W3", "netW1_sum", "netW2_sum", "netW3_sum")
# netW1_sum, netW2_sum, netW3_sum are summary measures...
# Define these using expressions, s.a., 
# W2[[0]] just means W2
# W2[[1:Kmax]] means vectors of (W2_netF_j), j=1, ..., Kmax
# sum(W2[[1:Kmax]]) is netW2_sum
# sW <- c(W1, W2, W3, sum(W1[ , 1:Kmax]), sum(W2[ , 1:Kmax]), sum(W3[ , 1:Kmax]))

# For deault, no need to define def.sW.gstar, only when requested, otherwise def.sW.gstar = def_sW_g0
def_sW_g0 <- def.sW(W1[[0]], W2[[0:Kmax]], W3[[0:Kmax]],
							netW1_sum = rowSums(W1[[1:Kmax]]),
							netW2_sum = rowSums(W2[[1:Kmax]]),
							netW3_sum = rowSums(W3[[1:Kmax]]))

def_sW_gstar <- def.sW(W1[[0]], W2[[0]], W3[[0]], type.g0 = FALSE)
# def_sW_gstar <- def.sW(W1 = W1, W2 = W2, W3 = W3, type.g0 = FALSE) # vs 2
# def_sW_gstar <- def.sW("W1[[0]]", "W2[[0]]", "W3[[0]]", type.g0 = FALSE) # vs 3
def_sA <- def.sA("A[[0:Kmax]]", sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))
# def_sA <- def.sA(A[[0:Kmax]], sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))

# NOT WRITTEN YET. Will evaluate the summary measures applied to the (O)bserved data (data.frame):
# res <- eval.summaries(summaries = def.summW.g0, Odata = df_K6, Kmax = kmax, NETIDnode = "Net_str", IDnode = "IDs")
# ....
# NOT WRITTEN YET. Need some kind of object that will hold intermediate inputs until tmlenet function if called...
# ....
# Qform <- "Y ~ I(netW2_1*(1-netA_1)+netW2_2*(1-netA_2)+netW2_3*(1-netA_3)+netW2_4*(1-netA_4)+netW2_5*(1-netA_5)+netW2_6*(1-netA_6))+I(netW3_1+netW3_2+netW3_3+netW3_4+netW3_5+netW3_6)"
(Qform <- "Y ~ I(" %+% paste(str_c(netvar("W2", (1:6)), "*", "(1-",netvar("A", (1:6)), ")"), collapse = "+") %+% ")" %+% "+"  %+% "I("  %+% paste(netvar("W3", (1:6)), collapse = "+") %+% ")")
# (Qform <- "Y ~ sum_1mAW2_nets + netW3_sum")
(gform <- "A ~  W1 + netW1_sum + netW2_sum + netW3_sum + nFriends")
(hform <- "sA ~ " %+% paste(netvar("W2", (1:6)), collapse = "+") %+% " + netW3_sum + nFriends")

#----------------------------------------------------------------------------------
# Example 1. Mean population outcome under deterministic intervention A=0 with 6 friends
# CAN RETIRE Qform, gform & hform -> fully replaced by summary measures
#----------------------------------------------------------------------------------
system.time(
tmlenet_K6out2 <- tmlenet(data = df_K6, Anode = "A", Wnodes = Wnodes, Ynode = "Y", nFnode = "nFriends",
						Kmax = kmax, IDnode = "IDs", NETIDnode = "Net_str", Qform = Qform, gform = gform, hform = hform,
						f_g1star = f.A_0, args_f_g1star = NULL, nQ.MCsims = 10, ng.MCsims = 10, alphaTMLE_B_only = FALSE,
						sW.g0 = def_sW_g0, sW.gstar = def_sW_gstar, sA = def_sA))
						# alternative way to pass summary measures:
						# sW.g0 = list("W1[[0]]", "W2[[0:Kmax]]", "W3[[0:Kmax]]", netW1_sum = "rowSums(W1[[1:Kmax]]"), netW2_sum = "rowSums(W2[[1:Kmax]])", netW3_sum = "rowSums(W3[[1:Kmax]])"), 
						# sW.gstar = NULL, 
						# sA = list("A[[0:Kmax]]", sum_1mAW2_nets = "rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))")
						
# N=1,000
# new h method has no subsetting (include degenerate outcome sA):
						   # old:     # new:
# epsilon (covariate)      0.02549743 0.04468718
# alpha (intercept)        0.05410938 0.09888994
# iptw epsilon (covariate) 0.03556655 0.03556655
			# old:		# new:
# tmle_A     0.5053725 0.5116754
# tmle_B     0.5051903 0.5119461
# iid.tmle_B 0.4475714 0.4475714
# tmle_iptw  0.5123310 0.5123310
# iptw_h     0.5065960 0.5185523
# iptw       0.4910014 0.4910014
# iid.iptw   0.4429414 0.4429414
# mle        0.4970377 0.4970377
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

# set.seed(543)
# df_K6 <-gendata_pop(nC=1, n_arr=2000, k_arr=kmax, EC_arr=EC, f.g_list="f.A", f.g_args_list=list(NULL))
# N=2,000
							# old:		# new:
# epsilon (covariate)      -0.01458956 -0.008449511
# alpha (intercept)        -0.03196439 -0.019919747
# iptw epsilon (covariate) -0.01916955 -0.019169554
			# old:		# new:
# tmle_A     0.4931451 0.4950080
# tmle_B     0.4931283 0.4948037
# iid.tmle_B 0.4379087 0.4379087
# tmle_iptw  0.4900913 0.4900913
# iptw_h     0.4816697 0.4865285
# iptw       0.4634407 0.4634407
# iid.iptw   0.4389123 0.4389123
# mle        0.4975765 0.4975765

# Verifying the old way of spec. Qform and new summary measures give the same answer:
# [1] "m.Q.init"
# Coefficients:
# (Intercept)  
# -0.9049
# I(W2_netF1 * (1 - A_netF1) + W2_netF2 * (1 - A_netF2) + W2_netF3 * (1 - A_netF3) + W2_netF4 * (1 - A_netF4) + W2_netF5 * (1 - A_netF5) + W2_netF6 * (1 - A_netF6))  
# 1.8026  
# I(W3_netF1 + W3_netF2 + W3_netF3 + W3_netF4 + W3_netF5 + W3_netF6)  
# -1.3216  
# Coefficients:
#    (Intercept)  sum_1mAW2_nets       netW3_sum  
#        -0.9049          1.8026         -1.3216  

tmlenet_K6out2$estimates$EY_g1.star$tmle_B
tmlenet_K6out2$estimates$EY_g1.star$CI_tmle_B_iidIC

#----------------------------------------------------------------------------------
# Example 2. Mean population outcome under deterministic intervention A=1 with 6 friends
#----------------------------------------------------------------------------------
tmlenet_K6out2 <- tmlenet(data=df_K6, Anode="A", Wnodes=Wnodes, Ynode="Y", nFnode="nFriends",
						Kmax=kmax, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform, h_form=hform,
						f.g1.star=f.A_1, f.g1_args=NULL, n_MCsims=10, n_samp_g0gstar=10)

tmlenet_K6out2$estimates$EY_g1.star$tmle_B
tmlenet_K6out2$estimates$EY_g1.star$CI_tmle_B_iidIC


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
#----------------------------------------------------------------------------------
tmlenet_out1 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", nFnode="nFriends",
						Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform, 
						f.g1.star=f.A_0, f.g1_args=NULL)

# TMLE estimate and iid IC-based 95% CI:
tmlenet_out1$estimates$EY_g1.star$tmle_B
tmlenet_out1$estimates$EY_g1.star$CI_tmle_B_iidIC

# Efficient IPTW (h) and iid IC-based 95% CI:
tmlenet_out1$estimates$EY_g1.star$iptw_h
tmlenet_out1$estimates$EY_g1.star$CI_iptw_h_iidIC

# Inefficient IPTW (g) + (two CIs, less conservative and more conservative)
tmlenet_out1$estimates$EY_g1.star$iptw
tmlenet_out1$estimates$EY_g1.star$CI_iptw_iidIC_1stO
tmlenet_out1$estimates$EY_g1.star$CI_iptw_iidIC_2ndO

# MLE
tmlenet_out1$estimates$EY_g1.star$mle

#----------------------------------------------------------------------------------
# Example 2. Mean population outcome under stochastic intervention P(A=1)=0.2
#----------------------------------------------------------------------------------
tmlenet_out2 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", nFnode="nFriends",
						Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
						f.g1.star=f.A_x, f.g1_args=list(x=0.2),
						n_MCsims=4000, n_samp_g0gstar=100)

# TMLE estimate and iid IC-based 95% CI:
tmlenet_out2$estimates$EY_g1.star$tmle_B
tmlenet_out2$estimates$EY_g1.star$CI_tmle_B_iidIC

# Efficient IPTW (h) and iid IC-based 95% CI:
tmlenet_out2$estimates$EY_g1.star$iptw_h
tmlenet_out2$estimates$EY_g1.star$CI_iptw_h_iidIC

# Inefficient IPTW (g) + (two CIs, less conservative and more conservative)
tmlenet_out2$estimates$EY_g1.star$iptw
tmlenet_out2$estimates$EY_g1.star$CI_iptw_iidIC_1stO
tmlenet_out2$estimates$EY_g1.star$CI_iptw_iidIC_2ndO

# MLE
tmlenet_out2$estimates$EY_g1.star$mle

#----------------------------------------------------------------------------------
# Example 3. Average treatment effect (ATE) for two interventions, f.g1.star: A=1 vs f.g2.star: A=0
#----------------------------------------------------------------------------------
tmlenet_out3 <- tmlenet(data=sample_network_k2, Anode="A", Wnodes="W1", Ynode="Y", nFnode="nFriends",
						Kmax=2, IDnode="IDs", NETIDnode="Net_str", Qform=Qform, gform=gform,
						f.g1.star=f.A_1, f.g1_args=NULL, f.g2.star=f.A_0, f.g2_args=NULL,
						n_MCsims=4000, n_samp_g0gstar=100)

# TMLE estimate for ATE + 95% CI
tmlenet_out3$estimates$ATE$tmle_B
tmlenet_out3$estimates$ATE$CI_tmle_B_iidIC

# Efficient IPTW (h) and iid IC-based 95% CI:
tmlenet_out3$estimates$ATE$iptw_h
tmlenet_out3$estimates$ATE$CI_iptw_h_iidIC

# Inefficient IPTW (g) + (two CIs, less conservative and more conservative)
tmlenet_out3$estimates$ATE$iptw
tmlenet_out3$estimates$ATE$CI_iptw_iidIC_1stO
tmlenet_out3$estimates$ATE$CI_iptw_iidIC_2ndO

# MLE
tmlenet_out3$estimates$ATE$mle



