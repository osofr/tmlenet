#***************************************************************************************
# LOAD DATA, LOAD A NETWORK
#***************************************************************************************
data(df_netKmax6) # load observed data
head(df_netKmax6)
data(NetInd_mat_Kmax6)  # load the network ID matrix
netind_cl <- simcausal:::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = 6)
netind_cl$NetInd <- NetInd_mat_Kmax6
head(netind_cl$nF)

#***************************************************************************************
# Example. Equivalent ways of defining the same summary measures.
# Note that 'nF' (number of friends per unit) is always added to def_sW() function call.
# Same rules apply to def_sA function, except that 'nF' is not added.
#***************************************************************************************
sW <- def_sW(W1, W2, W3)
sW <- def_sW("W1", "W2", "W3")
sW <- def_sW(W1 = W1, W2 = W2, W3 = W3)
sW <- def_sW(W1 = W1[[0]], W2 = W2[[0]], W3 = W3[[0]]) # W1[[0]] just means W1
sW <- def_sW(W1 = "W1[[0]]", W2 = "W2[[0]]", W3 = "W3[[0]]")

# evaluate the sW summary measures defined last:
resmatW <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
head(resmatW)

# define sA summary measures and evaluate:
sA <- def_sA(A, AW1 = A*W1)
resmatA <- sA$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
head(resmatA)

#***************************************************************************************
# Summary measures based on network (friend) values of the variable (matrix result).
#***************************************************************************************
# W2[[1:Kmax]] means vectors of W2 values of friends (W2_netF_j), j=1, ..., Kmax:
sW <- def_sW(netW2 = W2[[0:Kmax]], W3 = W3[[0]])
# evaluation result is a matrix:
resmat <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
# The mapping from the summary measure names to actual evaluation column names:
sW$sVar.names.map

# Equivalent way to define the same summary measure is to use syntax '+'
# and omit the names of the two summary measures above
# (the names are assigned automatically as "W2" for the first matrix W2[[0:Kmax]]
# and "W3" for the second summary measure "W3[[0]]")
sW <- def_sW(W2[[0:Kmax]]) + def_sW(W3[[0]])
resmat2 <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
head(resmat2)
# The mapping from the summary measure names to actual evaluation column names:
sW$sVar.names.map

#***************************************************************************************
# Define new summary measure as a sum of friend covariate values of W3:
#***************************************************************************************
# replaceNAw0 = TRUE sets all the missing values to 0
sW <- def_sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)

# evaluation result:
resmat <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)

#***************************************************************************************
# More complex summary measures that involve more than one covariate:
#***************************************************************************************
# replaceNAw0 = TRUE sets all the missing values to 0
sW <- def_sW(netW1W3 = W3[[1:Kmax]]*W3[[1:Kmax]])

# evaluation result (matrix):
resmat <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
# the mapping from the summary measure names to the matrix column names:
sW$sVar.names.map

#***************************************************************************************
# Vector results, complex summary measure (more than one unique variable name):
# NOTE: all complex summary measures must be named, otherwise an error is produced
#***************************************************************************************
# named expression:
sW <- def_sW(sum.netW2W3 = sum(W3[[1:Kmax]]*W2[[1:Kmax]]), replaceNAw0 = TRUE)
mat1a <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)

# the same unnamed expression (trying to run will result in error):
sW <- def_sW(sum(W3[[1:Kmax]]*W2[[1:Kmax]]), replaceNAw0 = TRUE)

\dontrun{
  mat1b <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)  
}

#***************************************************************************************
# Matrix result, complex summary measure (more than one unique variable name):
# NOTE: all complex summary measures must be named, otherwise an error is produced
#***************************************************************************************
# When more than one parent is present, the columns are named by convention:
# sVar.name%+%c(1:ncol)

# named expression:
sW <- def_sW(sum.netW2W3 = W3[[1:Kmax]]*W2[[1:Kmax]])
mat1a <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)

# the same unnamed expression (trying to run will result in error):
sW <- def_sW(W3[[1:Kmax]]*W2[[1:Kmax]])
\dontrun{
  mat1b <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
}

#***************************************************************************************
# Iteratively building higher dimensional summary measures using '+' function:
#***************************************************************************************
sW <- def_sW(W1) +
      def_sW(netW1 = W2[[1:Kmax]]) +
      def_sW(sum.netW1W3 = sum(W1[[1:Kmax]]*W3[[1:Kmax]]), replaceNAw0 = TRUE)

# resulting matrix of summary measures:
resmat <- sW$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
# the mapping from the summary measure names to the matrix column names:
sW$sVar.names.map

#***************************************************************************************
# Examples of summary measures defined by def_sA (functions of baseline and treatment)
#***************************************************************************************
sA <- def_sA(sum.netAW2net = sum((1-A[[1:Kmax]]) * W2[[1:Kmax]]),
            replaceNAw0 = TRUE) +
      def_sA(netA = A[[0:Kmax]])

resmat <- sA$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
sW$sVar.names.map

#***************************************************************************************
# More summary measures for sA
#***************************************************************************************
sA <- def_sA(netA = "A[[0:Kmax]]") +
      def_sA(sum.AW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0 = TRUE)

resmat <- sA$eval.nodeforms(data.df = df_netKmax6, netind_cl = netind_cl)
sW$sVar.names.map

#***************************************************************************************
# Using eval.summaries to evaluate summary measures for both, sW and sA
# based on the (O)bserved data (data.frame) and network
#***************************************************************************************
sW <- def_sW(netW2 = W2[[1:Kmax]]) +
      def_sW(netW3_sum = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)
            
sA <- def_sA(sum.AW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0 = TRUE) +
      def_sA(netA = A[[0:Kmax]])

data(df_netKmax6) # load observed data
data(NetInd_mat_Kmax6)  # load the network ID matrix
res <- eval.summaries(sW = sW, sA = sA, Kmax = 6, data = df_netKmax6,
                      NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

# Contents of the list returned by eval.summaries():
names(res)
# matrix of sW summary measures:
head(res$sW.matrix)
# matrix of sA summary measures:
head(res$sA.matrix)
# matrix of network IDs:
head(res$NETIDmat)
# Observed data (sW,sA) stored as "DatNet.sWsA" R6 class object:
res$DatNet.ObsP0
class(res$DatNet.ObsP0)
