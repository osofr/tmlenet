### --- Test setup ---
`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))


test.hforms_IPTW_positivity <- function() {
  data(datO_input)
  summary(datO_input)
  NetInd_mat <- attributes(datO_input)$netind_cl$NetInd_k
  K <- ncol(NetInd_mat)
  # tmlenet_options(useglm = FALSE)

  sW <-  def_sW(W1, W2, WNoise, corrW.F1, corrW.F2, corrW.F3, corrW.F4, corrW.F5, HUB = ifelse(nF >= 25, 1, 0), PA0 = (PA == 0))
  sA <-  def_sA(A, nF.PA = sum(PA[[1:Kmax]]), replaceNAw0 = TRUE) +
         def_sA(A.PAeq0 = A * (PA == 0)) +
         def_sA(nFPAeq0.PAeq1 = (nF.PA < 1) * (PA == 1)) +
         def_sA(sum.net.A = (sum(A[[1:Kmax]])*(HUB==0) + sum((W1[[1:Kmax]] > 4)*A[[1:Kmax]])*(HUB==1)),
                sum.net.A.sum.netPA = sum.net.A*nF.PA,
                replaceNAw0 = TRUE)

  new.sA1.stoch <-  def_new_sA(A = rbinom(n = length(A), size = 1, prob = 0.40))

  Qform <- "Y ~ nF.PA + A.PAeq0 + nFPAeq0.PAeq1 + sum.net.A + sum.net.A.sum.netPA + PA0 + W1 + W2 + corrW.F1 + corrW.F2 + corrW.F3 + corrW.F4 + corrW.F5"

  # This ordering gives bounded IPTW (was failing before):
  hform1.g0 <- "A + sum.net.A ~ HUB + nF.PA + nFPAeq0.PAeq1 + PA0"
  res1 <- tmlenet(data = datO_input, sW = sW, sA = sA,
                  Ynode = "Y", Kmax = K, NETIDmat = NetInd_mat,
                  intervene1.sA = new.sA1.stoch,
                  Qform = Qform, hform.g0 = hform1.g0, hform.gstar = hform1.g0,
                  optPars = list(
                    bootstrap.var = FALSE, n.bootstrap = 10
                    ))
  assertthat::assert_that(res1[["EY_gstar1"]][["estimates"]]["h_iptw",] < 0.20)

  # This ordering gives bounded IPTW, since A is now automatically put as a first covar when fitting sum.net.A:
  hform2.g0 <- "A + sum.net.A ~ HUB + nF.PA + PA0 + nFPAeq0.PAeq1"
  res2 <- tmlenet(data = datO_input, sW = sW, sA = sA,
                  Ynode = "Y", Kmax = K, NETIDmat = NetInd_mat,
                  intervene1.sA = new.sA1.stoch,
                  Qform = Qform, hform.g0 = hform2.g0, hform.gstar = hform2.g0,
                  optPars = list(
                    bootstrap.var = FALSE, n.bootstrap = 10
                    ))
  assertthat::assert_that(res2[["EY_gstar1"]][["estimates"]]["h_iptw",] < 0.20)

  # this ordering was working before (and works now as well), giving bounded IPTW:
  hform3.g0 <- "A + sum.net.A ~ HUB + PA0 + nF.PA + nFPAeq0.PAeq1"
  res3 <- tmlenet(data = datO_input, sW = sW, sA = sA,
                  Ynode = "Y", Kmax = K, NETIDmat = NetInd_mat,
                  intervene1.sA = new.sA1.stoch,
                  Qform = Qform, hform.g0 = hform3.g0, hform.gstar = hform3.g0,
                  optPars = list(
                    bootstrap.var = FALSE, n.bootstrap = 10
                    ))
  assertthat::assert_that(res3[["EY_gstar1"]][["estimates"]]["h_iptw",] < 0.20)
}