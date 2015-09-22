#***************************************************************************************
# Examples of summary measures defined by def.sW
#***************************************************************************************
# Use R expressions to define summary measure functions, s.a.,
# W2[[0]] just means W2; W2[[1:Kmax]] means vectors of (W2_netF_j), j=1, ..., Kmax;
def_sW <- def.sW(netW1 = W1[[0:Kmax]], noname = TRUE) +
          def.sW(netW2 = W2[[1:Kmax]], noname = TRUE) +
          def.sW(netW3 = W3[[1:Kmax]], noname = TRUE) +
          def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceNAw0 = TRUE)

def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
            def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceNAw0 = TRUE)

def_sW <- def.sW(netW2 = W2[[1:Kmax]], netW3_sum = rowSums(W3[[1:Kmax]]), noname= TRUE,
                  replaceNAw0 = TRUE)
def_sW <- def.sW(netW2 = W2[[1:Kmax]], netW3_sum = rowSums(W3[[1:Kmax]]),
                  replaceNAw0 = FALSE, noname = TRUE)

def_sW <- def.sW(W1 = W1[[0]], W1 = W2[[0:Kmax]], W1 = W3[[0:Kmax]],
                  netW1_sum = rowSums(W1[[1:Kmax]]),
                  netW2_sum = rowSums(W2[[1:Kmax]]),
                  netW3_sum = rowSums(W3[[1:Kmax]]))

#***************************************************************************************
# Examples of summary measures defined by def.sA
#***************************************************************************************
def_sA <- def.sA(sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]),
            replaceNAw0 = TRUE) +
            def.sA(netA = A[[0:Kmax]])

#***************************************************************************************
# Example. Three equivalent ways to define summary measures:
#***************************************************************************************
def_sW <- def.sW(W1 = W1[[0]], W2 = W2[[0]], W3 = W3[[0]])
def_sW <- def.sW(W1 = W1, W2 = W2, W3 = W3) # vs 2
def_sW <- def.sW(W1 = "W1[[0]]", W2 = "W2[[0]]", W3 = "W3[[0]]") # vs 3

#***************************************************************************************
# Example. See contents of the Define_sVar object returned by def.sW
#***************************************************************************************
def_sW$sVar.exprs
def_sW$sVar.expr.names
def_sW$ReplMisVal0
def_sW$sVar.misXreplace
def_sW$sVar.noname

#***************************************************************************************
# Example. More summary measures for sA
#***************************************************************************************
def_sA <- def.sA(netA = A[[0:Kmax]], noname = TRUE)
def_sA <- def.sA(netA = "A[[0:Kmax]]",
                sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))

#***************************************************************************************
# Example. More summary measures for sW and sA
#***************************************************************************************
def_sW <- def.sW(netW2 = W2[[1:Kmax]], noname = TRUE) + 
          def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceNAw0 = TRUE)
            
def_sA <- def.sA(sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]),
            replaceNAw0 = TRUE) +
          def.sA(netA = A[[0:Kmax]], noname = TRUE)

#***************************************************************************************
# Evaluate the summary measures applied to the  (O)bserved data (data.frame) and network
#***************************************************************************************
data(df_netKmax6) # load observed data
data(NetInd_mat_Kmax6)  # load the network ID matrix
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

