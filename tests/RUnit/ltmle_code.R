# install.packages("simcausal")
# ************************************************************************************
# tmlenet version of ltmle needs to be installed from github:
# devtools::install_github('osofr/tmlenet', ref = "ltmle", build_vignettes = FALSE)
# actual ltmle with 2 versions of pooling:
# devtools::install_github('joshuaschwab/ltmle', ref = "pooling", build_vignettes = FALSE)
# ************************************************************************************

if(FALSE) {
  library("RUnit")
  library("roxygen2")
  library("devtools")
  setwd(".."); setwd(".."); getwd()
  document()
  load_all("./", create = FALSE) # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # tmlenet:::debug_set() # SET TO DEBUG MODE
  # setwd("..");
  # install("tmlenet", build_vignettes = FALSE) # INSTALL W/ devtools:
}


# simcausal DAG for sim1:
make.simcausal.DAG <- function(t.end=16){
  require(simcausal)
  options(simcausal.verbose = FALSE)
  Drm <- DAG.empty()
###### L(0)=(Z(0),I*(0)) - I*(m) represents value of test at m-1 if test performed: I*(m)=I(m)*T(m-1) where I(m) is value of test (possibly unobserved)
  Drm <- Drm +
    node("Y", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  #Z(0)
    node("lastNat1", t = 0, distr = "rconst", const = 0) +  #Z(0) - see below for definition,  set at 0 at t = 0 by convention (see also below)
    node("highA1c.UN", t = 0, distr = "rbern", prob = 0.05) +  # I(0) - possibly unobserved lab
    node("highA1c", t = 0, distr = "rbern", prob = highA1c.UN[0]) +  # I*(0) = I(0)*T(-1) with T(-1) = 1 by convention - all patient come with an A1c value known - T(-1) is what I call N(-1)
    node("CVD", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 1, 0.5, 0.1)) +  #Z(0)
    node("timelowA1c.UN", t = 0, distr = "rbern", prob = 1-highA1c.UN[0]) +  # Z(0) - counts the number of time the (possibly unobserved) A1c was low
###### A(0)
    node("TI", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 0, ifelse(CVD[0] == 1, 0.5, 0.1), ifelse(CVD[0] == 1, 0.9, 0.5))) +
    # node("C", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  # no censoring in first bin
###### T(0)
    node("N", t = 0, distr = "rbern", prob = 1) +
###### L(t) for t = 1, ..., t.end
    # node("Y", t = 1:t.end, distr = "rbern", prob = plogis(-6.5 + 1*CVD[0] + 4*highA1c.UN[t-1] + 0.05*timelowA1c.UN[t-1]),  EFU = TRUE) +  # Z(t)
    node("Y", t = 1:t.end, distr = "rbern", prob = ifelse(Y[t-1]==1,1L,plogis(-6.5 + 1*CVD[0] + 4*highA1c.UN[t-1] + 0.05*timelowA1c.UN[t-1]))) +  # Z(t)
    node("lastNat1", t = 1:t.end, distr = "rconst", const = ifelse(N[t-1] == 0, lastNat1[t-1] + 1, 0)) +  # Z(1)  just a function of past \bar{N}(t-1) - 0 probs current N at 1,  1 probs previous N a 1,  2 probs  the one before the previous was at 1,  etc.
    node("highA1c.UN", t = 1:t.end, distr = "rbern", prob = ifelse(TI[t-1] == 1, 0.1, ifelse(highA1c.UN[t-1] == 1, 0.9, min(1, 0.1 + t/.(t.end))))) +  # I(t)
    node("highA1c", t = 1:t.end, distr = "rbern", prob = ifelse(N[t-1] == 1, highA1c.UN[t], highA1c[t-1])) + # I*(m)=I(m)*T(m-1)  (I actually replace I*(m)=0 with when T(m-1)=0 with last value carried forward,  i.e. I*(m-1)
    node("timelowA1c.UN", t=1:t.end, distr="rnorm", mean=sum(1-highA1c.UN[0:t]),  sd=0) +  # Z(m)
###### A(t)
    node("TI", t = 1:t.end, distr = "rbern",
      prob =
        ifelse(TI[t-1] == 1, 1,
          ifelse(N[t-1] == 1,
            ifelse(highA1c[t] == 0,
              ifelse(CVD[0] == 1, 0.3, 0.1),
                ifelse(CVD[0] == 1, 0.7, 0.5)), 0))) +
    # node("C", t = 1:t.end, distr = "rbern", prob = ifelse(t == t.end, 1, 0), EFU = TRUE) +  # last time point is when admin end of study occurs
###### N(t)
    node("N", t = 1:t.end, distr = "rbern", prob = 1)
  return(Drm)
}

# ------------------------------------------------------------------------------------
# big-data simcausal-based simulation for sequential G-COMPUTATION formula
# ------------------------------------------------------------------------------------
sim1.simcausal <- function() {
  t.end <- 50 # number of time-points
  n <- 100000 # sample size
  abar <- 0   # intervention

  # ----------------------------------------------------------------------------
  # Generate some data
  # ----------------------------------------------------------------------------
  require(simcausal)
  dagobj <- set.DAG(make.simcausal.DAG(t.end=t.end))
  simdat <- sim(dagobj, n = 100000)
  (allnodes <- names(attributes(simdat)$node_nms))
  simdatfit <- simdat[,-c(1,2)]
  rem.vars <- c("highA1c.UN_"%+%0:t.end, "timelowA1c.UN_"%+%0:t.end, "lastNat1_"%+%0:t.end, "highA1c_"%+%t.end,"TI_"%+%t.end,"N_"%+%t.end)
  simdatfit <- simdatfit[,!names(simdatfit)%in%rem.vars]
  names(simdatfit)
  mean(simdatfit[,"Y_"%+%t.end])

  # define nodes and Qform regressions:
  nodes <- CreateNodes(simdatfit, Anodes = "TI_"%+%1:(t.end-1), Cnodes = NULL,
                       Lnodes = c("highA1c_"%+%2:(t.end-1), "Y_"%+%2:(t.end-1)),
                       Ynodes="Y_"%+%t.end)
  # Qform <- CreateLYNodes(simdatfit, nodes, check.Qform=TRUE, Qform=Qform)$Qform
  Qform.default <- GetDefaultForm(simdatfit, nodes, is.Qform=TRUE)
  # lapply(Qform.default, RhsVars)

  # define intervention (matrix of counterfactual treatment assignments):
  Abar.df <- data.frame(matrix(abar, nrow=n, ncol=length(1:(t.end-1))))
  colnames(Abar.df) <- "TI_"%+%1:(t.end-1)
  dim(Abar.df)
  head(Abar.df)

  # Define input data objects for tmlenet:
  def_sW <- def.sW("CVD_0",
                    highA1c=highA1c[0:(t.end-1)],
                    N=N[0:(t.end-1)],
                    Y=Y[1:(t.end-1)])
  datnetW <- DatNet$new()$make.sVar(Odata = simdatfit, sVar.object = def_sW)

  def_sA <- def.sA(TI=TI[1:(t.end-1)])
  datnetA <- DatNet$new()$make.sVar(Odata = simdatfit, sVar.object = def_sA)
  # observed data object:
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = simdatfit[,nodes$Y])$make.dat.sWsA()
  # counterfactual abar data object:
  datnetA.gstar <- DatNet$new()$make.sVar(Odata = Abar.df, sVar.object = def_sA)
  datNet.gstar <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA.gstar)$make.dat.sWsA()
  # head(datNet.gstar$dat.sVar)

  # mean(simdatfit[,"Y_"%+%t.end])
  # t(head(cbind(simdatfit[,"Y_"%+%t.end],simdatfit[,nodes$Y]),100))

  # Set-up sequential G-comp regs:
  newreg_Gcomp <- SGCompRegClass$new(outvar = "Y", nodes = nodes, Qforms = Qform.default)
  newsummaries_Gcomp <- SGcompSummariesModel$new(reg = newreg_Gcomp)
  # fit the G-comp regression:
  timerun <- system.time(
    fitted_Gcomp <- newsummaries_Gcomp$fit(data = datNetObs, newdata = datNet.gstar)
  )
  print(timerun/60)
  #      user    system   elapsed
  # 7.8283000 0.6985333 8.5077667

  head(datNetObs$getQ.kplus1())
  mean(datNetObs$getQ.kplus1())
  mean(data[,"Y"])
}

# ------------------------------------------------------------------------------------
# ltmle simulation from Josh
# ------------------------------------------------------------------------------------
sim2.validation <- function() {
  #uses github branch: "pooling"
  # devtools::install_github('joshuaschwab/ltmle', ref = "pooling", build_vignettes = FALSE)
  library(ltmle)
  rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
  GenerateData <- function(n, est.psi0=FALSE) {
    W <- rbinom(n, 1, 0.5)
    if (est.psi0) {
      A1 <- rep(1, n)
    } else {
      A1 <- rexpit(-3 + W)
    }
    L <- rexpit(W + A1)
    if (est.psi0) {
      A2 <- rep(1, n)
    } else {
      A2 <- rexpit(-1.5 + W + A1 + L)
    }
    Y <- rexpit(W + k*A1 + L + A2)
    if (est.psi0) return(mean(Y))
    return(data.frame(W, A1, L, A2, Y))
  }

  k <- -4
  n <- 50000
  data <- GenerateData(n, FALSE)
  head(data)
  abar <- c(1, 1)

  # RunIter <- function(j) {
  j <- 1
  set.seed(j+100000)
  data <- GenerateData(n, FALSE)
  est.strat <- rep(NA, 2)
  est.nostrat <- rep(NA, 2)
  names(est.strat) <- names(est.nostrat) <- c("pooling1 (set last Abar)", "set all Abars")
  
  # Qform <- c(L="Q.kplus1~W*A1", Y="Q.kplus1~W*A1*L*A2") #saturated models, gives means within strata, estimates for method1 and method2 are identical
  Qform <- c(L="Q.kplus1~W+A1", Y="Q.kplus1~W+A1+L+A2") #different in finite samples

  # ----------------------------------------------------------------------------------------
  # Running 4 GCOMP versions in ltmle package by pooling method and stratification on abar
  # ----------------------------------------------------------------------------------------
  for (i in 1:2) {
    pooling1 <<- i == 1
    gform <- matrix(1, nrow(data), 2) #not used in gcomp, but saves a little time in iptw
    timerun.ltmle1 <- system.time(
      r.strat <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y",
                  abar=abar, gcomp=TRUE, IC.variance.only = TRUE,
                  estimate.time = FALSE, gform=gform, stratify = TRUE, Qform=Qform)
    )
    est.strat[i] <- r.strat$estimates[1]
    # print(timerun.ltmle1)
    timerun.ltmle2 <- system.time(
      r.nostrat <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y",
                  abar=abar, gcomp=TRUE, IC.variance.only = TRUE,
                  estimate.time = FALSE, gform=gform, stratify = FALSE, Qform=Qform)
    )
    est.nostrat[i] <- r.nostrat$estimates[1]
  }
  print(timerun.ltmle2)
  #  user  system elapsed
  # 2.674   0.315   2.976

  # print(timerun.ltmle2)
  print("estimates for stratify=TRUE"); print(est.strat)
  print("estimates for stratify=FALSE"); print(est.nostrat)
  # cat(est, est[1]-est[2], "\n")
  # return(est)
# }
# for (j in 1:20) RunIter(j)

  # ----------------------------------------------------------------------------------------
  # Running GCOMP from tmlenet that corresponds to ltmle with stratify=FALSE and pooling1=FALSE
  # (pool all observations when fitting Q, set all Abars for the prediction step)
  # ----------------------------------------------------------------------------------------
  nodes <- CreateNodes(data, Anodes = c("A1", "A2"), Cnodes = NULL, Lnodes = c("L"), Ynodes="Y")
  Abar.df <- data.frame(matrix(abar, nrow=n, ncol=2, byrow=TRUE))
  colnames(Abar.df) <- c("A1", "A2")

  # Define data:
  def_sW <- def.sW("W", "L")
  def_sA <- def.sA("A1", "A2")
  datnetW <- DatNet$new()$make.sVar(Odata = data, sVar.object = def_sW)
  datnetA <- DatNet$new()$make.sVar(Odata = data, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = data[,nodes$Y])$make.dat.sWsA()
  # Counterfactual trt assignments (for prediction):
  datnetA.gstar <- DatNet$new()$make.sVar(Odata = Abar.df, sVar.object = def_sA)
  datNet.gstar <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA.gstar)$make.dat.sWsA()
  # head(datNet.gstar$dat.sVar)

  # Set-up sequential G-comp regs:
  newreg_Gcomp <- SGCompRegClass$new(outvar = "Y", nodes = nodes, Qforms = Qform)
  newsummaries_Gcomp <- SGcompSummariesModel$new(reg = newreg_Gcomp)
  # Run G-comp estimation:
  timerun <- system.time(
    fitted_Gcomp <- newsummaries_Gcomp$fit(data = datNetObs, newdata = datNet.gstar)
  )
  print(timerun)
  #  user  system elapsed
  # 0.119   0.041   0.160
  timerun.ltmle2/timerun
  #      user    system   elapsed
  # 22.470588  7.682927 18.600000

  gcomp.psi.n <- mean(datNetObs$getQ.kplus1())
  print("tmlenet gcomp.psi.n: " %+% gcomp.psi.n)
  cat("difference between ltmle & tmlenet G-comps: ", gcomp.psi.n-est.nostrat[2])
}

# ------------------------------------------------------------------------------------
# modified ltmle simulation from Josh adding another time-point
# ------------------------------------------------------------------------------------
sim3.validation <- function() {
  #uses github branch: "pooling"
  # devtools::install_github('joshuaschwab/ltmle', ref = "pooling", build_vignettes = FALSE)
  library(ltmle)
  rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
  GenerateData <- function(n, est.psi0=FALSE) {
    W <- rbinom(n, 1, 0.5)
    if (est.psi0) {
      A1 <- rep(1, n)
    } else {
      A1 <- rexpit(-3 + W)
    }
    L1 <- rexpit(W + A1)
    if (est.psi0) {
      A2 <- rep(1, n)
    } else {
      A2 <- rexpit(-1.5 + W + A1 + L1)
    }
    L2 <- rexpit(W + A1 + L1 + A2)
    if (est.psi0) {
      A3 <- rep(1, n)
    } else {
      A3 <- rexpit(-1.5 + W + A1 + A2 + 0.4*L1 + 0.4*L2)
    }
    L3 <- rexpit(W + A1 + A2 + 0.4*A3 + L1 + 0.4*L2)
    if (est.psi0) {
      A4 <- rep(1, n)
    } else {
      A4 <- rexpit(-1.5 + W + A1 + A2 + A3 + L1 + L2 + L3)
    }
    Y <- rexpit(W + L1 + 0.3*L2 + 0.2*L3 + k*A1 + 0.4*A2 + 0.4*A3 + 0.1*A4)
    if (est.psi0) return(mean(Y))
    return(data.frame(W, A1, L1, A2, L2, A3, L3, A4, Y))
  }

  k <- -4
  n <- 500000
  data <- GenerateData(n, FALSE)
  head(data)
  abar <- c(1, 1, 1, 1)

  # RunIter <- function(j) {
  j <- 1
  set.seed(j+100000)
  data <- GenerateData(n, FALSE)
  est.strat <- rep(NA, 2)
  est.nostrat <- rep(NA, 2)
  names(est.strat) <- names(est.nostrat) <- c("pooling1 (set last Abar)", "set all Abars")
  
  # Qform <- c(L="Q.kplus1~W*A1", Y="Q.kplus1~W*A1*L*A2") #saturated models, gives means within strata, estimates for method1 and method2 are identical
  Qform <- c(L1="Q.kplus1~W+A1",
            L2="Q.kplus1~W+A1+L1+A2",
            L3="Q.kplus1~W+A1+A2+A3+L1+L2",
            Y="Q.kplus1~W+A1+A2+A3+A4+L1+L2+L3")

  # ----------------------------------------------------------------------------------------
  # Running 4 GCOMP versions in ltmle package by pooling method and stratification on abar
  # ----------------------------------------------------------------------------------------
  # for (i in 1:2) {
  i <- 2
    pooling1 <<- i == 1
    gform <- matrix(1, nrow(data), ncol=4) #not used in gcomp, but saves a little time in iptw
    timerun.ltmle1 <- system.time(
      r.strat <- ltmle(data, Anodes=c("A1","A2","A3","A4"), Lnodes=c("L1","L2","L3"), Ynodes="Y",
                  abar=abar, gcomp=TRUE, IC.variance.only = TRUE,
                  estimate.time = FALSE, gform=gform, stratify = TRUE, Qform=Qform)
    )
    est.strat[i] <- r.strat$estimates[1]
    # print(timerun.ltmle1)
    timerun.ltmle2 <- system.time(
      r.nostrat <- ltmle(data, Anodes=c("A1","A2","A3","A4"), Lnodes=c("L1","L2","L3"), Ynodes="Y",
                  abar=abar, gcomp=TRUE, IC.variance.only = TRUE,
                  estimate.time = FALSE, gform=gform, stratify = FALSE, Qform=Qform)
    )
    est.nostrat[i] <- r.nostrat$estimates[1]
  # }
  print(timerun.ltmle2)
   #   user  system elapsed
   # 40.590   5.281  45.805
  print("estimates for stratify=TRUE"); print(est.strat)
  print("estimates for stratify=FALSE"); print(est.nostrat)
  # return(est)
  # }
  # for (j in 1:20) RunIter(j)

  # ----------------------------------------------------------------------------------------
  # Running GCOMP from tmlenet that corresponds to ltmle with stratify=FALSE and pooling1=FALSE
  # (pool all observations when fitting Q, set all Abars for the prediction step)
  # ----------------------------------------------------------------------------------------
  nodes <- CreateNodes(data, Anodes = c("A1","A2","A3","A4"), Cnodes = NULL, Lnodes = c("L1","L2","L3"), Ynodes="Y")
  Abar.df <- data.frame(matrix(abar, nrow=n, ncol=4, byrow=TRUE))
  colnames(Abar.df) <- c("A1","A2","A3","A4")

  # Define data:
  def_sW <- def.sW("W", "L1","L2","L3")
  def_sA <- def.sA("A1", "A2", "A3", "A4")
  datnetW <- DatNet$new()$make.sVar(Odata = data, sVar.object = def_sW)
  datnetA <- DatNet$new()$make.sVar(Odata = data, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = data[,nodes$Y])$make.dat.sWsA()
  # Counterfactual trt assignments (for prediction):
  datnetA.gstar <- DatNet$new()$make.sVar(Odata = Abar.df, sVar.object = def_sA)
  datNet.gstar <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA.gstar)$make.dat.sWsA()
  # head(datNet.gstar$dat.sVar)

  # Set-up sequential G-comp regs:
  newreg_Gcomp <- SGCompRegClass$new(outvar = "Y", nodes = nodes, Qforms = Qform)
  newsummaries_Gcomp <- SGcompSummariesModel$new(reg = newreg_Gcomp)
  # Run G-comp estimation:
  timerun <- system.time(
    fitted_Gcomp <- newsummaries_Gcomp$fit(data = datNetObs, newdata = datNet.gstar)
  )
  print(timerun)
  #  user  system elapsed
  # 3.256   0.730   3.974
  timerun.ltmle2/timerun
  #      user    system   elapsed
  # 14.520502  5.087266 11.730072
  gcomp.psi.n <- mean(datNetObs$getQ.kplus1())
  print("tmlenet gcomp.psi.n: " %+% gcomp.psi.n)
  cat("difference between ltmle & tmlenet G-comps: ", gcomp.psi.n-est.nostrat[2])

}



# Organize nodes
CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes) {
  Anodes <- NodeToIndex(data, Anodes)
  Cnodes <- NodeToIndex(data, Cnodes)
  Lnodes <- NodeToIndex(data, Lnodes)
  Ynodes <- NodeToIndex(data, Ynodes)
  nodes <- list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, AC=sort(c(Anodes, Cnodes)))
  nodes$LY <- CreateLYNodes(data, nodes, check.Qform=FALSE)
  return(nodes)
}
# Convert named nodes to indicies of nodes
NodeToIndex <- function(data, node) {
  if (! is.data.frame(data)) stop("data must be a data frame")
  if (is.numeric(node) || is.null(node)) return(node)
  if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
  index <- match(node, names(data))
  if (any(is.na(index))) {
    stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
  }
  return(index)
}
# Get the LY nodes but don't include "blocks" of L/Y nodes uninterrupted by A/C nodes
CreateLYNodes <- function(data, nodes, check.Qform, Qform) {
  LYnodes <- sort(c(nodes$L, nodes$Y))
  #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
  nodes.to.remove <- NULL
  if (length(LYnodes) > 1) {
    for (i in 1:(length(LYnodes) - 1)) {
      cur.node <- LYnodes[i]
      next.node <- LYnodes[i + 1]
      if (! any(cur.node:next.node %in% nodes$AC)) {
        nodes.to.remove <- c(nodes.to.remove, next.node)
      }
    }
  }
  new.LYnodes <- setdiff(LYnodes, nodes.to.remove)
  if (check.Qform) {
    removed.Qform.index <- NULL
    for (i in nodes.to.remove) {
      index <- which(names(Qform) == names(data)[i])
      if (length(index) > 0) {
        removed.Qform.index <- c(removed.Qform.index, index)
      }
    }
    if (! is.null(removed.Qform.index)) {
      message("L/Y nodes (after removing blocks)  : ", names(data)[new.LYnodes], "\n")
      message("Qform names                        : ", names(Qform), "\n")
      message(paste("The following nodes are not being considered as L/Y nodes because they are part of a block of L/Y nodes. They are being dropped from Qform:\n"), paste(names(Qform)[removed.Qform.index], "\n", collapse=" "))
      Qform <- Qform[-removed.Qform.index]
    }
    return(list(LYnodes=new.LYnodes, Qform=Qform))
  }
  return(new.LYnodes)
}
# Get the default Q or g formula - each formula consists of all parent nodes except censoring and event nodes [also except A nodes if stratifying]
GetDefaultForm <- function(data, nodes, is.Qform, stratify=FALSE, survivalOutcome=TRUE, showMessage=TRUE) {
  if (is.Qform) {
    lhs <- rep("Q.kplus1", length(nodes$LY))
    node.set <- nodes$LY
  } else {
    lhs <- names(data)[nodes$AC]
    node.set <- nodes$AC
  }
  if (stratify) {
    stratify.nodes <- c(nodes$C, nodes$A)
  } else {
    stratify.nodes <- c(nodes$C)
  }
  if (survivalOutcome) {
    stratify.nodes <- c(stratify.nodes, nodes$Y)
  }
  form <- NULL
  for (i in seq_along(node.set)) {
    cur.node <- node.set[i]
    if (cur.node == 1) {
      form[i] <- paste(lhs[i], "~ 1")  #no parent nodes
    } else {
      parent.node.names <- names(data)[setdiff(1:(cur.node - 1), stratify.nodes)]
      if (length(parent.node.names) == 0) {
        form[i] <- paste(lhs[i], "~ 1")
      } else {
        form[i] <- paste(lhs[i], "~", paste(parent.node.names, collapse=" + "))
      }
    }
    names(form)[i] <- names(data)[cur.node]
  }
  
  if (showMessage) {
    #Prints formulas with automatic wrapping thanks to print.formula
    message(ifelse(is.Qform, "Qform", "gform"), " not specified, using defaults:")
    lapply(seq_along(form), function(i, names) {
      message("formula for ", names[i], ":")
      #Using print on a formula because it nicely wraps
      message(capture.output(print(as.formula(form[i]), showEnv=FALSE)))
    }, names=names(form))
    message("")
  }
  return(form)
}
