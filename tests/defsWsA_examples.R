#----------------------------------------------------------------------------------
# EXAMPLES OF SUMMARY MEASURES
#----------------------------------------------------------------------------------
# Defined with R expressions, s.a., 
# W2[[0]] just means W2; W2[[1:Kmax]] means vectors of (W2_netF_j), j=1, ..., Kmax; sum(W2[[1:Kmax]]) is netW2_sum
def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
            def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceMisVal0 = TRUE)
def_sA <- def.sA(sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]), replaceMisVal0 = TRUE) +
            def.sA(netA = A[[0:Kmax]])

# MORE EXAMPLES of various summary measures:
def_sW <- def.sW(netW1 = W1[[0:Kmax]], noname = TRUE) +
          def.sW(netW2 = W2[[1:Kmax]], noname = TRUE) +
          def.sW(netW3 = W3[[1:Kmax]], noname = TRUE) +
          def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceMisVal0 = TRUE)

def_sW <- def.sW(W2[[1:Kmax]], netW3_sum = rowSums(W3[[1:Kmax]]), replaceMisVal0 = TRUE)
def_sW <- def.sW(W2[[1:Kmax]], netW3_sum = rowSums(W3[[1:Kmax]]), replaceMisVal0 = FALSE)
def_sW <- def.sW(W1[[0]], W2[[0:Kmax]], W3[[0:Kmax]],
                  netW1_sum = rowSums(W1[[1:Kmax]]),
                  netW2_sum = rowSums(W2[[1:Kmax]]),
                  netW3_sum = rowSums(W3[[1:Kmax]]))

#----------------------------------------------------------------------------------
# Example. 3 equivalent ways to define summary measures:
#----------------------------------------------------------------------------------
def_sW <- def.sW(W1 = W1[[0]], W2 = W2[[0]], W3 = W3[[0]])
def_sW <- def.sW(W1 = W1, W2 = W2, W3 = W3) # vs 2
def_sW <- def.sW(W1 = "W1[[0]]", W2 = "W2[[0]]", W3 = "W3[[0]]") # vs 3

#----------------------------------------------------------------------------------
# Example. See contents of the Define_sVar object
#----------------------------------------------------------------------------------
def_sW$sVar.exprs
def_sW$sVar.expr.names
def_sW$ReplMisVal0
def_sW$sVar.misXreplace
def_sW$sVar.noname

#----------------------------------------------------------------------------------
# Example. More summary measures for sA
#----------------------------------------------------------------------------------
def_sA <- def.sA(netA = A[[0:Kmax]], noname = TRUE)
def_sA <- def.sA(netA = "A[[0:Kmax]]", sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]))

#----------------------------------------------------------------------------------
# Example. More summary measures for sW and sA
#----------------------------------------------------------------------------------
def_sW <- def.sW(netW2 = W2[[1:Kmax]], noname = TRUE) +
            def.sW(netW3_sum = rowSums(W3[[1:Kmax]]), replaceMisVal0 = TRUE)
            
def_sA <- def.sA(sum_1mAW2_nets = rowSums((1-A[[1:Kmax]]) * W2[[1:Kmax]]), replaceMisVal0 = TRUE) +
            def.sA(netA = A[[0:Kmax]], noname = TRUE)

#----------------------------------------------------------------------------------
# Example.  NOT WRITTEN YET. Will evaluate the summary measures applied to the (O)bserved data (data.frame):
# res <- eval.summaries(summaries = def_sA, Odata = df_K6, Kmax = kmax, NETIDnode = "Net_str", IDnode = "IDs")
#----------------------------------------------------------------------------------