#***************************************************************************************
# Define some summary measures for sW and sA
#***************************************************************************************
def_sW <- def.sW(netW2 = W2[[1:Kmax]], noname = TRUE) + 
          def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceNAw0 = TRUE)
            
def_sA <- def.sA(sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]),
            replaceNAw0 = TRUE) +
          def.sA(netA = A[[0:Kmax]], noname = TRUE)

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

