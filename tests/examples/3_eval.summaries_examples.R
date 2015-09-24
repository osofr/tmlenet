#***************************************************************************************
# Define some summary measures for sW
#***************************************************************************************
def_sW <- def.sW(W1, W2, W3) +
          def.sW(netW1 = W1[[1:Kmax]]) +
          def.sW(mean.netW2 = mean(W2[[1:Kmax]]), replaceNAw0 = TRUE) +
          def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)

#***************************************************************************************
# Define some summary measures for sA
#***************************************************************************************
def_sA <- def.sA("A") + # equivalent to def.sA(A)
          def.sA(netA = A[[0:Kmax]]) +
          def.sA(sum.AW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0 = TRUE)

#***************************************************************************************
# Evaluate the summary measures applied to the  (O)bserved data (data.frame) and network
#***************************************************************************************
# load observed data:
data(df_netKmax6)
# load the network ID matrix:
data(NetInd_mat_Kmax6)
res <- eval.summaries(sW = def_sW, sA = def_sA,  Kmax = 6, data = df_netKmax6,
  NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

#***************************************************************************************
# Contents of the list returned by eval.summaries():
#***************************************************************************************
names(res)
# matrix of sW summary measures:
res$sW.matrix
# matrix of sA summary measures:
res$sA.matrix
# matrix of network IDs:
res$NETIDmat
# Observed data (sW,sA) stored as "DatNet.sWsA" R6 class object:
res$DatNet.ObsP0
class(res$DatNet.ObsP0)

