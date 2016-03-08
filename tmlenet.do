
#------------------------------------------------------------------------------------------------
# tmlenet TO DO LIST
#------------------------------------------------------------------------------------------------

- Pre-save all simulated test datasets. Current method of setting seed and running simcausal is unstable.
- Put parserfunctions_R6.r back into simcausal (Setting opts$NoChangeFunCalls <- FALSE in simcausal)
- Export Define_sVar from simcausal, import simcausal::Define_sVar into tmlenet namespace
- See what other global vars in parserfunctions_R6.r will need to be redefined in tmlenet

# Standardize the naming conventions across package:
	- Remove "." from all function names (replace with _)
	- Remove Upper case from all function names

-tmlenet: Implement data-adaptive weight truncation wrt minimization of MSE (taking max of 5% weight of total as truth)

-tmlenet: Allow SL to fit Q_N, g_N and h (i.e P(A_j|A_{1},..,A_{j-1}, W_i\\inF_j))

-tmlenet: Pooled P(sA|sW) estimation:
	Need to validate pooled P(sA|sW) estimation for continuous sA with in iid case! 
	Need example where pooled model is true -> use IPTW for estimating g, pooled over time (use RN sim)

-tmlenet 17: +5 "R/tmlenet.R" * 00/00/00 00:00
	Add option to set fixed seed for MC sampling (so that two tmlenet runs with MC sim w/ same inputs and same seed return the same estimates)

-tests 17: +10 "tests/RUnit_tests_03_iidcont_sA_tests.R" * 00/00/00 00:00
	verify iid estimation (no network, no nF) still works in tests when calling tmlenet()

-tmlenet, interface 17: +5 "R/tmlenet.R" * 00/00/00 00:00
	Specifying exposure model for static f_gstar: for most-nonparametric example of sA should not be fitting glms, 
	since P(sA|sW) under f_gstar is just indicators

-tmlenet, interface 17: +5 "R/tmlenet.R" * 00/00/00 00:00
	Re-create optPars list at the top of the tmlenet() call and print message with current tmlenet settings

-tmlenet, m.Q.init 17: +5 "R/tmlenet.R" * 15/05/23 00:00
	In the future add option to fit separate m.Q.init model for each value of nF?

-tmlenet, m.Q.init 17: +5 "R/tmlenet.R" * 15/05/23 00:00
	Move fitting of m.Q.init inside get_all_ests?

-fit.hbars 17: +5 "R/tmlenet.R" * 15/05/23 00:00
	Currently enforcing the outcome summary in h.g0.sVars$outvars and h.gstar.sVars$outvars to be the same. In the future might allow different summary measures for sA_nms_g0 and sA_nms_gstar.

-tmlenet 17: +5 "R/tmlenet.R" * 15/06/29 23:48
	pre-evaluate summary measures on small batch of data to get dims of sA & sW and to check for errors

-DatNet, DatNet.sWsA 28: +0 "R/DatNetClass.R" * 15/07/01 13:26
	is nFnode even needed DatNet anymore???? Where is it used? Can it be removed completely?

-DatNet.sWsA 28: +0 "R/DatNetClass.R" * 15/07/01 13:26
	rename datnetW, datnetA to O.datnetW, O.datnetA for clarity
	*** NOTE *** When sVar is cat might be better to set bin_bymass = FALSE to avoid collapsing of categories for sVar
	Fix 71 will still not solve the issue for ordinal sVar. Need to manually set intervals when sVar is categorical!

-DatNet, DatNet.sWsA 18: +0 "R/DatNet_class.R" * 15/07/01 14:10
	(OPTIONAL) ENABLE ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO DatNet
	Need to save det node flags as a separate mat, can't add them to sVar since all sVars will be automatically added to A ~ predictors	

-get_all_ests 19: +0 "R/tmlenet.R" * 15/06/30 20:51
	use glm.fit or speedglm.Wfit for m.Q.star, m.Q.star_reg_A, m.Q.star_reg_B, m.Q.star_iptw

-get_all_ests 22: +0 "R/tmlenet.R" * 15/06/30 20:53
	rename get_all_ests into fit.param.fluct

-get_all_ests 23: +0 "R/tmlenet.R" * 15/06/30 20:53
	move MC estimation (and h estimation?) outside get_all_ests

-ContinSummaryModel 23: +0 "R/SummariesModelClass.R" * 15/06/30 20:53
	(low priority) see how to generalize to pooled fits, k-specific fits, etc (use subset definitions + ?)

-ContinSummaryModel 23: +0 "R/SummariesModelClass.R" * 15/06/30 20:53
	(***BUG***) Currently make.bins_mtx_1 fails on binary A with automatic bin detection (all values are placed in last (2nd) bin)

-CategorSummaryModel 23: +0 "R/SummariesModelClass.R" * 15/06/30 20:53
	Need to test that make.bins_mtx_1 will do the right thing when x_cat is (0, ..., ncats) instead of (1, ..., ncats)

-CategorSummaryModel 23: +0 "R/SummariesModelClass.R" * 15/06/30 20:53
	Add to categorical exposure tests a call to tmlenet (not just iptw)

-SummariesModel 23: +0 "R/SummariesModelClass.R" * 15/06/30 20:53
	rename all predict() to predictP_1
	rename all predictAeqa() to predict.like.P_a

-DefineSummariesClass 23: +0 "R/DefineSummariesClass.R" * 15/06/30 20:53

-sVar_evaluator, sVar.name 23: +0 "R/DefineSummariesClass.R" * 15/06/30 20:53
	Validate the test that checks the resulting column names in sVar are always unique!

-sVar_evaluator 28: +0 "R/sVar_evaluator.R" * 15/07/01 13:26
	Consider returning sVar.res_l instead of mat.sVar, also see if there are faster alternatives to cbind (i.e., pre-allocating sVar.mat); perform benchmarks to see if there is any noticable benefit

-sVar_evaluator 28: +0 "R/sVar_evaluator.R" * 15/07/01 13:26

-+..DefineSummariesClass 28: +0 "R/sVar_evaluator.R"  * 15/07/01 13:26
  Allow adding character vector summary measures for sVar2, s.a., 
  def_sW(W2[[1:Kmax]]) + "netW3_sum = rowSums(W3[[1:Kmax]]"
  Accept sA & sW as character vectors / lists passed to tmlenet (in addition to current set-up)
  When sW / sA are just lists of character vectors need to capture the calling env and call DefineSummariesClass constructor:
    # user.env <- parent.frame()
    # user.env_l <- list(user.env = user.env)
    # sW <- do.call(DefineSummariesClass$new, c(sW, list(type = "sW"), user.env_l))
    # sW.gstar <- do.call(DefineSummariesClass$new, c(sW.gstar, list(type = "sW.gstar"), user.env_l))
    # sA <- do.call(DefineSummariesClass$new, c(sA, list(type = "sA"), user.env_l))  

-BinDat, BinOutModel 28: +0 "R/DatNetClass.R" * 15/07/01 13:26
	(Low priority) Consider merging these two classes (BinDat,BinOutModel) into one
