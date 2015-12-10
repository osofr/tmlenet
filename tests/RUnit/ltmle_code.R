# install.packages("simcausal")
# ************************************************************************************
# ltmle needs to be installed from github:
# devtools::install_github('osofr/tmlenet', ref = "ltmle", build_vignettes = FALSE)
# ************************************************************************************


# if(FALSE) {
#   library("RUnit")
#   library("roxygen2")
#   library("devtools")
#   setwd(".."); setwd(".."); getwd()
#   document()
#   load_all("./", create = FALSE) # load all R files in /R and datasets in /data. Ignores NAMESPACE:
#   # tmlenet:::debug_set() # SET TO DEBUG MODE
#   # setwd("..");
#   # install("tmlenet", build_vignettes = FALSE) # INSTALL W/ devtools:
# }


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
# simcausal-based simulation for sequential G-COMPUTATION formula
# ------------------------------------------------------------------------------------
sim1.simcausal <- function() {

  t.end <- 50 # number of time-points
  dagobj <- set.DAG(make.simcausal.DAG(t.end=t.end))
  simdat <- sim(dagobj, n = 100000)
  (allnodes <- names(attributes(simdat)$node_nms))
  simdatfit <- simdat[,-c(1,2)]
  rem.vars <- c("highA1c.UN_"%+%0:t.end, "timelowA1c.UN_"%+%0:t.end, "lastNat1_"%+%0:t.end, "highA1c_"%+%t.end,"TI_"%+%t.end,"N_"%+%t.end)
  simdatfit <- simdatfit[,!names(simdatfit)%in%rem.vars]
  names(simdatfit)
  mean(simdatfit[,"Y_"%+%t.end])

  nodes <- CreateNodes(simdatfit, Anodes = "TI_"%+%1:(t.end-1), Cnodes = NULL,
                       Lnodes = c("highA1c_"%+%2:(t.end-1), "Y_"%+%2:(t.end-1)),
                       Ynodes="Y_"%+%t.end)
  # Qform <- CreateLYNodes(simdatfit, nodes, check.Qform=TRUE, Qform=NULL)$Qform
  Qform.default <- GetDefaultForm(simdatfit, nodes, is.Qform=TRUE)
  lapply(Qform.default, RhsVars)

  # Define data for tmlenet package:
  def_sW <- def.sW("CVD_0",
                    highA1c=highA1c[0:(t.end-1)],
                    N=N[0:(t.end-1)],
                    Y=Y[1:(t.end-1)])
  datnetW <- DatNet$new()$make.sVar(Odata = simdatfit, sVar.object = def_sW)
  def_sA <- def.sA(TI=TI[1:(t.end-1)])
  datnetA <- DatNet$new()$make.sVar(Odata = simdatfit, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = simdatfit[,nodes$Y])$make.dat.sWsA()

  mean(simdatfit[,"Y_"%+%t.end])
  t(head(cbind(simdatfit[,"Y_"%+%t.end],simdatfit[,nodes$Y]),100))

  # Set-up sequential G-comp regs:
  newreg_Gcomp <- SGCompRegClass$new(outvar = "Y", nodes = nodes, Qforms = Qform.default)
  newsummaries_Gcomp <- SGcompSummariesModel$new(reg = newreg_Gcomp)
  # fit the G-comp regression:
  timerun <- system.time(
    fitted_Gcomp <- newsummaries_Gcomp$fit(data = datNetObs)
  )
  print(timerun/60)

  head(datNetObs$getQ.kplus1())
  mean(datNetObs$getQ.kplus1())
  mean(data[,"Y"])
}



# ------------------------------------------------------------------------------------
# ltmle simulation from Josh
# ------------------------------------------------------------------------------------
sim2 <- function() {
  #uses github branch: "pooling"
  rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
  GenerateData <- function(n, est.psi0=FALSE) {
    W <- rbinom(n, 1, 0.5)
    if (est.psi0) {
      A_1 <- rep(1, n)
    } else {
      A_1 <- rexpit(-3 + W)
    }
    L1_1 <- rexpit(W + A_1)
    L2_1 <- rexpit(W + A_1)
    if (est.psi0) {
      A_2 <- rep(1, n)
    } else {
      A_2 <- rexpit(-1.5 + W + A_1 + L1_1)
    }
    Y <- rexpit(W + k*A_1 + L1_1 + A_2)
    if (est.psi0) return(mean(Y))
    return(data.frame(W, A_1, L1_1, L2_1, A_2, Y))
  }
  k <- -4
  n <- 1000
  data <- GenerateData(n, FALSE)

  names(data)
  nodes <- CreateNodes(data, Anodes = c("A_1", "A_2"), Cnodes = NULL, Lnodes = c("L1_1", "L2_1"), Ynodes="Y")
  Qform <- CreateLYNodes(data, nodes, check.Qform=TRUE, Qform=NULL)$Qform
  Qform.default <- GetDefaultForm(data, nodes, is.Qform=TRUE)
  lapply(Qform.default, RhsVars)

  # Define data:
  def_sW <- def.sW("W", "L1_1", "L2_1")
  def_sA <- def.sA("A_1", "A_2")
  datnetW <- DatNet$new()$make.sVar(Odata = data, sVar.object = def_sW)
  datnetA <- DatNet$new()$make.sVar(Odata = data, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = data[,nodes$Y])$make.dat.sWsA()

  # Set-up sequential G-comp regs:
  newreg_Gcomp <- SGCompRegClass$new(outvar = "Y", nodes = nodes, Qforms = Qform.default)
  newsummaries_Gcomp <- SGcompSummariesModel$new(reg = newreg_Gcomp)
  fitted_Gcomp <- newsummaries_Gcomp$fit(data = datNetObs)
  head(datNetObs$getQ.kplus1())
  mean(datNetObs$getQ.kplus1())
  mean(data[,"Y"])

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

# # ltmleMSM for a single final.Ynode
# FixedTimeTMLE <- function(inputs) {
#   inputs$summary.measures <- NULL #just to make sure it isn't used - should only use combined.summary.measures
#   data <- inputs$data
#   nodes <- inputs$nodes
#   num.regimes <- 1

#   fit.Qstar <- vector("list", length(nodes$LY))
#   names(fit.Qstar) <- names(data)[nodes$LY]
#   Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow=n, ncol=num.regimes)
#   for (LYnode.index in length(nodes$LY):1) {
#     cur.node <- nodes$LY[LYnode.index]
#     # subsetting by intervention, will be a separate function:
#     deterministic.list.origdata <- IsDeterministic(data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
#     uncensored <- IsUncensored(data, nodes$C, cur.node)
#     intervention.match <- subs <- matrix(nrow=n, ncol=num.regimes)
#     i <- 1
#     abar <- GetABar(inputs$regimes, i)
#     intervention.match[, i] <- InterventionMatch(data, abar=abar, nodes$A, cur.node)
#     newdata <- SetA(data, abar=abar, nodes, cur.node)
#     deterministic.list.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
#     if (inputs$stratify) {
#       subs[, i] <- uncensored & intervention.match[, i] & !deterministic.list.origdata$is.deterministic
#     } else {
#       subs[, i] <- uncensored & !deterministic.list.origdata$is.deterministic
#     }

#     if (any(subs[, i])) {
#       Q.est <- Estimate(inputs$Qform[LYnode.index],
#                         data=data.frame(data,
#                         Q.kplus1=Qstar.kplus1[, i]),
#                         family="quasibinomial",
#                         newdata=newdata,
#                         subs=subs[, i],
#                         type="link",
#                         nodes=nodes)
#       logitQ[, i] <- Q.est$predicted.values
#     } else {
#       if (! all(deterministic.list.newdata$is.deterministic)) {
#         msg <- paste0("ltmle failed trying to estimate ", inputs$Qform[LYnode.index], " because there are no observations that are\nuncensored", ifelse(inputs$stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.Q.function\n")
#         stop(msg)
#       }
#       Q.est <- list(fit="no estimation of Q at this node because all rows are set deterministically")
#     }
#   }
#   return(list(IC=IC, Qstar=Qstar, cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, g.ratio=g.ratio, est.var=est.var, fit=list(g=fit.g, Q=fit.Q, Qstar=fit.Qstar))) 
# }





