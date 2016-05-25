#***************************************************************************************
# Define some summary measures for sW
#***************************************************************************************
sW <- def_sW(W1, W2, W3) +
      def_sW(netW1 = W1[[1:Kmax]]) +
      def_sW(netW2 = W2[[1:Kmax]]) +
      def_sW(mean.netW2 = mean(W2[[1:Kmax]]), replaceNAw0 = TRUE) +
      def_sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)

#***************************************************************************************
# Define some summary measures for sA
#***************************************************************************************
sA <- def_sA(A = A, netA = A[[0:Kmax]]) +
      def_sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0 = TRUE)

#***************************************************************************************
# Evaluate the summary measures applied to the  (O)bserved data (data.frame) and network
#***************************************************************************************
# load observed data:
data(df_netKmax6)
# load the network ID matrix:
data(NetInd_mat_Kmax6)
res <- eval.summaries(sW = sW, sA = sA,  Kmax = 6, data = df_netKmax6,
  NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

#***************************************************************************************
# Contents of the list returned by eval.summaries():
#***************************************************************************************
names(res)
# observed data matrix with (sW,sA) summary measures:
head(res$DatNet.ObsP0$dat.sWsA)
# matrix of network IDs:
head(res$NETIDmat)
# Observed data summary measures (sW,sA) and network
# stored as "DatNet.sWsA" R6 class object:
res$DatNet.ObsP0
class(res$DatNet.ObsP0)

#***************************************************************************************
# Using DatNet.ObsP0 as input to tmlenet():
#***************************************************************************************
options(tmlenet.verbose = FALSE) # set to TRUE to print status messages

res_K6 <- tmlenet(data =  df_netKmax6, NETIDmat = NetInd_mat_Kmax6,
                    Kmax = 6,
                    sA = sA, sW = sW,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    Ynode = "Y",
                    intervene1.sA = def_new_sA(A = 1L))

# res_K6 <- tmlenet(DatNet.ObsP0 = res$DatNet.ObsP0,
#                     Qform = "Y ~ sum.netW3 + sum.netAW2",
#                     hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
#                     hform.gstar = "netA ~ sum.netW3",
#                     Ynode = "Y",
#                     intervene1.sA = def_new_sA(A = 0L))
# f_gstar1 = 0L)
# Anodes = "A",

res_K6$EY_gstar1$estimates
res_K6$EY_gstar1$IC.vars
res_K6$EY_gstar1$IC.CIs



