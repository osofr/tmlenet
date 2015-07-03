# Uncomment and configure. ALL matched lines matters:
# file [absolute_path/]filename.ext
# mysql 127.0.0.1 username password scheme
# http typetodo.com repository [username password]

-simfcts 1: +0 "../../../../Dropbox/KP/monitoring_simstudy/sim1_using_simcausal/simfcts.r" **Anon 15/07/02 11:36
	TO ADD: if (noNtint), create a model based on true TI.t that ignores N.t

-sVar_evaluator, sVar.name 2: +0 "R/sVar_evaluator.R" **Anon 15/06/29 14:10
	Instead of throwing an error for un-named vector results (i.e., def.sW.g0(A)), would be nice to extract the first is.name of the expression and use that as a name.

-sVar_evaluator, sVar.name 3: +0 "R/sVar_evaluator.R" **Anon 15/06/29 14:10
	Allow overwriting names for some ncol(expr_res) > 1 in sVar.res_l. Allow using naturally assigned names for some ncol(expr_res) > 1 in sVar.res_l

-sVar_evaluator, sVar.name 4: +0 "R/sVar_evaluator.R" **Anon 15/06/28 11:50
	Create an option $keep.sVar.nms; When TRUE do not change the output column names for sVar mat with ncol > 1. Create an option to overwrite sVar mat colname(s) with user-provided names

-sVar_evaluator, UI 5: +1 "R/sVar_evaluator.R" **Anon 15/06/28 11:51
	Need a diagnostic tool that will evaluate and return the result of the summary measures applied to user Odata data.frame...

-sVar_evaluator, NetInd_k 6: +1 "R/sVar_evaluator.R" **Anon 15/06/29 14:09
	TEST NetInd_k, generating sA under g_star sampled data.df: nrow(NetInd_k) = n, while nrow(data.df) = n*p  => COULD BE A PROBLEM. NEED TO CHECK WORKS AS INTENDED.

-sVar_evaluator, sVar.name 7: +0 "R/sVar_evaluator.R" **Anon 15/06/29 14:11
	Check that the resulting column names in sVar are all unique!

-tmlenet 8: +0 "" **Anon 15/06/29 21:02
	check all names exist in data (Anode, Wnodes, Ynode, etc...)

+tmlenet 9: +0 "R/tmlenet.R" **Anon 15/06/29 16:08
	if no gform / hform specified, create a default (all names in sW.g0 + nF?)

-tmlenet 10: +0 "" **Anon 15/06/29 21:02
	need default setting for sW.g0, sW.gstar & sA when these are left unscpecified

-tmlenet 11: +0 "" **Anon 15/06/29 21:02
	Wnodes & Anode are no longer needed when sW, sA, give a warning that Wnodes,Anodes will be ignored

+DatNetClass 12: +0 "R/DatNet_class.R" **Anon 15/07/01 14:21
	**** THIS IS DANGEROUS AND NEEDS TO BE CHANGED **** Removing self$df.sVar <- df.sW.sA since at the moment df.sVar is an active binding in DatNet.

-BinOutModel, predict, predictAeqa 13: +0 "tests/runittests.R" **Anon 15/06/30 10:53
	Need to be linked together, since can create a discrepancy for missing(newdata) but !missing(obs.DatNet.sWsA)

+DatNet.sWsA 14: +0 "R/tmlenet.R" **Anon 15/07/01 09:13
	DatNet.sWsA field IS TO BE RENAMED TO O.datnetA for clarity ***

-BinOutModel, predictAeqa 15: +0 "tests/runittests.R" **Anon 15/06/30 10:53
	TO DO: MOVE PART OF THE CODE TO self

-tmletnet 16: +0 "R/tmlenet.R" **Anon 15/06/29 23:48
	make below part of DatNet evaluation

-tmletnet 17: +5 "R/tmlenet.R" **Anon 15/06/29 23:48
	pre-evaluate summary measures on small batch of data to get dims of sA & sW and to check for errors

-DatNet, DatNet.sWsA 18: +0 "R/DatNet_class.R" **Anon 15/07/01 14:10
	(OPTIONAL) ENABLE ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO DatNet

-get_all_ests 19: +0 "R/tmlenet.R" **Anon 15/06/30 20:51
	use glm.fit or speedglm.Wfit for m.Q.star_reg_A

-get_all_ests 20: +0 "R/tmlenet.R" **Anon 15/06/30 20:52
	use glm.fit or speedglm.Wfit for m.Q.star_reg_B

-get_all_ests 21: +0 "R/tmlenet.R" **Anon 15/06/30 20:52
	use glm.fit or speedglm.Wfit for m.Q.star_iptw

-get_all_ests 22: +0 "R/tmlenet.R" **Anon 15/06/30 20:53
	rename get_all_ests into fit.param.fluct

-get_all_ests 23: +0 "R/tmlenet.R" **Anon 15/06/30 20:53
	move MC estimation (and h estimation?) outside get_all_ests

-tmlenet 24: +0 "R/tmlenet.R" **Anon 15/07/01 08:39
	need to sort out how to store sW.gstar and sW.g0 and how to pass those to fit.hbars

-DatNet, add_det 25: +0 "R/DatNet_class.R" **Anon 15/07/01 14:10
	Need to save det node flags as a separate mat, can't add them to sVar since all sVars will be automatically added to A ~ predictors

-ContinSummaryModel$new 26: +0 "R/Summaries_class.R" **Anon 15/07/01 12:00
	consider launching O.datnetA$set.sVar.intrvls(...) when intervals for sVar aren't defined

-ContinSummaryModel 27: +0 "R/Summaries_class.R" **Anon 15/07/01 12:01
	Put subset eval in a separate function (with var as arg + additional args) +

-sVar_evaluator 28: +0 "R/sVar_evaluator.R" **Anon 15/07/01 13:26
	Consider returning sVar.res_l instead of mat.sVar, also see if there are faster alternatives to cbind (i.e., pre-allocating sVar.mat); perform benchmarks to see if there is any noticable benefit

-sVar_evaluator / misXreplace 29: +0 "R/DatNet_class.R" **Anon 15/07/01 14:10
	Might have to set-up sVar-specific misValRepl = TRUE/FALSE in make.sVar

+DatNet.sWsA 30: +0 "R/DatNet_class.R" **Anon 15/07/01 14:45
	NEW 07/01/15: For datnetAstar.sVar to work it needs to have access to Odata where A is replaced with A^*

-evalsubst 31: +0 "R/DatNet_class.R" **Anon 15/07/01 15:44
	NEED TO RE-WRITE ALL SUBSET EVALUATION SO THAT IT WORKS WITH MATRICES (can't use expressions for subset anymore)

-simfcts 32: +0 "../../../../Dropbox/KP/monitoring_simstudy/sim1_using_simcausal/simfcts.r" **Anon 15/07/02 11:36
	fit a propenisty score model for true P(A.t=1|L.t) that ignores N(t-1)

