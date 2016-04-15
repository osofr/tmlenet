# --------------------------------------------------------------------------------------------
# NON-INFORMATIVE ERROR FOR PARAMETRIC BOOTSTRAP
# for a node in boot.nodes = c("PA","A") not having a regression specified in either in hform.g0 or boot.regs
# results in error: "Error: attempt to apply non-function"
# need to catch that and make more informative
# --------------------------------------------------------------------------------------------

library("tmlenet")
library("simcausal")
options(tmlenet.verbose = FALSE)

# --------------------------------------------------------------------------------------------
# Specificy the observed data distribution with a DAG object
# --------------------------------------------------------------------------------------------
create_OnetDAG <- function() {
  `%+%` <- function(a, b) paste0(a, b)
  require("igraph")
  require("simcausal")
  options(simcausal.verbose=FALSE)
  # --------------------------------------------------------------------------------------------
  # preferential attachment (BA) model (power law deg distribution):
  # --------------------------------------------------------------------------------------------
  generate.igraph.prefattach <- function(n, power, zero.appeal, m, ...) {
    g <- sample_pa(n, power = power, zero.appeal = zero.appeal, m = m)
    g <- as.directed(as.undirected(g, mode = "collapse"), mode = "mutual")
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(g)
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    return(NetInd_out$NetInd_k)
  }
  # --------------------------------------------------------------------------------------------
  # Defining the DAG object
  # --------------------------------------------------------------------------------------------
  D <- DAG.empty()
  D <- D +
      network("Net.prefattach", netfun = "generate.igraph.prefattach", power = 0.5, zero.appeal = 5, m = 5) +
      node("LatentW", distr = "rcat.b0", probs = c(0.0494, 0.1823, 0.2806, 0.2680,0.1651, 0.0546)) +
      node("HUB", distr = "rconst", const = ifelse(nF >= 25, 1, 0)) # is this person a hub?
  D <- D +
      node("nF", distr = "rconst", const = nF) + # number of friends
      node("W1", distr = "rcat.b1", probs = c(0.0494, 0.1823, 0.2806, 0.2680,0.1651, 0.0546)) +
      node("W2", distr = "rbern", prob = plogis(-0.2)) +
      node("WNoise", distr = "rbern", prob = plogis(-0.4)) +
      # FOR IID INDEPENDENT W's:
      node("corrW.F1", distr = "rbern", prob = plogis(-8 + 2*LatentW + 2*LatentW)) +
      node("corrW.F2", distr = "rbern", prob = plogis(-6 + 1.5*LatentW + 1.5*LatentW)) +
      node("corrW.F3", distr = "rbern", prob = plogis(-6 + 1.5*LatentW + 1.5*LatentW)) +
      node("corrW.F4", distr = "rbern", prob = plogis(-4 + LatentW + LatentW)) +
      node("corrW.F5", distr = "rbern", prob = plogis(-4 + LatentW + LatentW))
  D <- D +
      node("PA", distr = "rbern", prob = W2*0.05 + (1-W2)*0.15) + # Physically active at baseline (depends on W2)
      node("nF.PA", distr = "rconst", const = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) # number of phys. active friends
  D <- D + node("A", distr = "rbern", prob = 0.25)
  D <- D + node("sum.net.A", distr = "rconst", const = (sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1)), replaceNAw0 = TRUE)
  D <- D +
      node("probY", distr = "rconst",
          const = plogis(ifelse(PA==1,
                  +5 - 15*(nF.PA < 1), # the probability of retaining PA drops if no friends are PA
                  -8.0 + 0.25*A) +
                  +0.5*sum.net.A + 0.25*nF.PA*sum.net.A + 0.5*nF.PA +
                  +0.5*(W1-1) - 0.58*W2 +
                  -0.5*(3.477-1) + 0.58*0.4496 +
                  +4*corrW.F1 -2*corrW.F2 -2*corrW.F3 +2*corrW.F4 -2*corrW.F5 +
                  -4*0.6841   +2*0.6727   +2*0.6724   -2*0.6513   +2*0.6538),
          replaceNAw0 = TRUE)

  D <- D + node("Y", distr = "rbern", prob = probY)
  D <- set.DAG(D, latent.v = c("LatentW"), n.test = 200)
}


Dset <- create_OnetDAG()
net.seed <- 123345 # net.seed <- NULL
rndseed.reset.node <- "LatentW" # to condition on the network (same network is sampled each time)
datO <- sim(Dset, n = 1000, rndseed = net.seed, rndseed.reset.node = rndseed.reset.node)
netind_cl <- attributes(datO)$netind_cl
NetInd_mat <- attributes(datO)$netind_cl$NetInd
K <- ncol(NetInd_mat)

hform.g0 <- hform.gstar <- "nF.PA ~ HUB + nF"
Qform <- "Y ~ nF.PA + A.PAeq0 + nFPAeq0.PAeq1 + sum.net.A + sum.net.A.sum.netPA + PA0 + W1 + W2 + corrW.F1 + corrW.F2 + corrW.F3 + corrW.F4 + corrW.F5"

sW <-  def_sW(W1, W2, WNoise, corrW.F1, corrW.F2, corrW.F3, corrW.F4, corrW.F5, HUB = ifelse(nF >= 25, 1, 0), PA0 = (PA == 0))
sA <-  def_sA(A, nF.PA = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) +
       def_sA(A.PAeq0 = A * (PA == 0)) +
       def_sA(nFPAeq0.PAeq1 = (nF.PA < 1) * (PA == 1)) +
       def_sA(sum.net.A = (sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1)),
              sum.net.A.sum.netPA = sum.net.A*nF.PA,
              replaceNAw0 = TRUE)

new.sA1.nFPA_A.8b <- def_new_sA(nF.PA = ifelse((nF <= 15) | (W1 <= 2), nF.PA + 1, nF.PA)) +
                     def_new_sA(A = rbinom(n = length(A), size = 1, prob = ifelse(nF >= 20, 0.9, ifelse(nF >= 15, 0.40, 0))))

new.sA2.nFPA_A <- def_new_sA(nF.PA = nF.PA) +
                  def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.1))


checkException(
  )
res <- tmlenet(data = datO, sW = sW, sA = sA,
              Ynode = "Y",
              Kmax = K, NETIDmat = NetInd_mat,
              intervene1.sA = new.sA1.nFPA_A.8b,
              intervene2.sA = new.sA2.nFPA_A,
              Qform = Qform,
              hform.g0 = hform.g0,
              hform.gstar = hform.gstar,
              optPars = list(
                bootstrap.var = TRUE,
                n.bootstrap = 10,
                boot.nodes = c("PA","A"),
                boot.form = c("PA ~ W2")
                )
              )

