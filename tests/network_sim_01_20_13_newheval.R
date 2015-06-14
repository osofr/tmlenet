rm(list=ls())
setwd("/Users/monikaizano/Dropbox/F'12_Research/Simulation")
getwd()

library(locfit)
library(xtable)
library(bigmemory)
library(biganalytics)
library(parallel)
library(plyr)
# library(data.table)
# library(Matrix)
# library(tmle)
#source("tmle_ISM.R")
options(bigmemory.typecast.warning=FALSE)

# simulate # of network friends from uniform (vs fixed)
NONFIXED_N_FRIENDS <- T

# processors <- detectCores()
processors <- 1
# setwd("./")
print("No. of CPUs")
print(processors)

#-------------------------------------------------------------------------------------------------------
  # Define coefficient for i individual's W
#-------------------------------------------------------------------------------------------------------
  # coefficients for true g(A|W), for k=5
  Intercept_A <- -0.5 
  coeff_A_W <- 2.5
  coeff_A_Wfriend <- 0.9 
  # coefficients for true Q(Y|A,W), k=5
  Intercept_Y <- -2.5 
  coeff_Y_W <- 2
  coeff_Y_Wfriend <- 0.7 
  coeff_Y_A <- 3
  coeff_Y_Afriend <- 1.0 
  
#-------------------------------------------------------------------------------------------------------
# Community Level Interventions g*C(A|W)
#------------------------------------------------------------------------------------------------------- 
# Actual pop g(A|W) from NPSEM, logistic fcn with a coeff for network Ws and a coeff for individ W_i
	  	# DEFINITION OF g(A_i|W)=P[A=1|cA]:      
  f.A <- function(W1, cA_Pa, n, ...) {				
	  	f.netWs <- function(netwk_W, coeff)  sum(coeff * netwk_W)     	#sum all coeffs for network W
	  	sum_friendWs <- sapply(cA_Pa, f.netWs, coeff_A_Wfriend)
	  	indivW <- coeff_A_W * W1      									#individual's covariate W_i
	  	# print("Parameters for g(A|W):")
	  	# print(paste("Intercept_A= ", Intercept_A, "coeff_A_W= ", coeff_A_W))
	  	# F(x) = 1 / (1 + exp(-(x-m)/s))
	  	pA.man <- 1 / (1+exp(-(Intercept_A + indivW + sum_friendWs)))
	 	return(pA.man)  }
# Deterministically setting every A=0
  f.A_0 <- function(n, ...) rep(0,n)		
# Deterministically setting every A=1
  f.A_1 <- function(n, ...) rep(1,n)		
# Deterministically set A=1 based on cutt-off value for Ws or if W_i = 1
  f.A_cutt_offW <- function(W1, cA_Pa, n, cutt_offW, ...) {	
	  	f.netWs <- function(netwk_W)  sum(netwk_W) 	
	  	sum_friendWs <- sapply(cA_Pa, f.netWs)
	  	indivW <- W1
	  	A <- rep(0,n)
	  	A[which((indivW==1)|(sum_friendWs>cutt_offW))] <- 1 
	 	return(A)  	}
# Determin to A=1 based on connectivity |N_i|= {low, high}, based on median
  f.A_1_highNi <- function(n, nFriends, ...) {		
	    Ni_med = quantile(nFriends, 0.5)
	    #print(paste("Ni_low", Ni_med))
	  	A <- rep(0,n)
		A[which(nFriends<= Ni_med)] <- 0
		A[which(nFriends> Ni_med)] <- 1
	 	return(A)  	 } 	
# Only personal covariates, no interference - DIRECT EFFECT, g(A|W_i)
  f.A_nospill <- function(W1, n, ...) plogis(Intercept_A + coeff_A_W*W1)
# Only interference effects A_i - INDIRECT EFFECT, g(A|W\W_i)
  f.A_spillonly <- function(W1, cA_Pa, n, ...) {	
	  	f.netWs <- function(netwk_W, coeff)  sum(coeff * netwk_W)
	  	sum_friendWs <- sapply(cA_Pa, f.netWs, coeff_A_Wfriend)
	  	pA <- plogis(Intercept_A + sum_friendWs)
	 	return(pA)  	  }
# Set x% of community to A=1
  f.A_x <- function(n, x, ...) rep(x, n) 
# Set x% of community to A=1 stratified by W_i
  f.A_xlevelW <- function(W1, n, x_Wi_0, x_Wi_1, ...) {		
	  	x <- rep(0,n)
		x[which(W1==0)] <- x_Wi_0
		x[which(W1==1)] <- x_Wi_1
	  	return(x)   } 	 
# Set x% of community to A=1 based on connectivity |N_i|= {low, high}, based on median
  f.A_xlevelNi <- function(n, nFriends, x_Ni_low, x_Ni_high, ...) {		
	    Ni_med = quantile(nFriends, 0.5)
	    #print(paste("Ni_low", Ni_med))
	  	x <- rep(0,n)
		x[which(nFriends<= Ni_med)] <- x_Ni_low
		x[which(nFriends> Ni_med)] <- x_Ni_high
	  	return(x)  } 	
 
 #-------------------------------------------------------------------------------------------------------
# Generate the network population using structural equations model
#-------------------------------------------------------------------------------------------------------
#Sample 1 community (C_j) with EC and f.g_A_W=A^C
gendata_Cj <- function(C_j = 1, n, k, EC, f.g_name=NULL, f.g_args=NULL) { 
  	#n - number of individuals
  	#k - max size of each individual's network
  	#C_j - community # being sampled (j)
  	#EC - community-level covariate, only influences W_i (right now its a probability for binom distr of W_1)
  	#f.g_name - fcn for community intervention on A^C, i.e. g(A|W)
  	#f.g_args - additional args to be passed to f.g_name (as a list)
  	#----------------------------------------------------------------
  	# Defining structural equations 
  	#----------------------------------------------------------------
  	.f.W1 <- function(Cj_prob, n) rbinom(n, 1, Cj_prob)
  	.f.Net_num <- function(samp=FALSE, n, k) if (samp) sample(0:k, n, replace=T) else rep(k,n)	
  	# Sample connectivity matrix (0,1), 
  	# randomly selecting from the list of available individuals 
  	# so that we never exceed |N_i| for each i 
  	# (the algorithm may result in actual of friends being < |N_i| for some individuals)
  	.f.genConnectMatx <- function(Net_num) {	
	  	I <- big.matrix(n,n, type="short", init=0)
		nFriendTot <- array(rep(0,n))
	  	for (index in (1:n)) {
		  		I[index,index] <- 1
		  		FriendSampSet <- setdiff(which(nFriendTot<Net_num), c(0:index)) #set of possible friends to sample	
				# if ((index%%1000)==0) print(paste("index",index))		#output every 1,000th index count
				nFriendSamp <- max(Net_num[index] -	nFriendTot[index], 0)	#check i's network is not already filled to max
				nFriendSamp <- min(nFriendSamp, length(FriendSampSet)) #make sure i still has someone to connect with 
				if ((length(FriendSampSet)==1) & (nFriendSamp>0))  { 
					friends_i <- FriendSampSet	#to fix the situation when FriendSampSet=n_i - one number
					} else {
						friends_i <- sort(sample(FriendSampSet, nFriendSamp))}	#sample from the possible friend set
				if (nFriendSamp>0) {
					I[index, friends_i] <- 1 
					I[friends_i, index] <- 1
					nFriendTot[index] <- nFriendTot[index] + nFriendSamp
					nFriendTot[friends_i] <- nFriendTot[friends_i] + 1 
				}
		}
		return(I)	
	}
  	#Update the # of friends for each individual (given sampled connnectivity matrix)
  	.f.Net_num_update <- function(ConnectMatx) return(colsum(ConnectMatx, c(1:n)) - 1)	
  	#Convert connectivity matx to a vector of friend's ids for each individual i, packaged together in a list
  	.f.Net_vect <- function(ConnectMatx) {	  	
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(which(ConnectMatx[index,]!=0), index)
				netwklist_i	}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}
  	.f.Net_vect_big <- function(ConnectMatx) {#Same, but using mwhich() instead of which() (faster for n>30,000)
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
				netwklist_i	}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}		
  	.f.mkstrNet <- function(Net) sapply(Net, function(Net_i) paste(Net_i, collapse=" "))	
  	.f_g_wrapper <- function(fcn_name, ...) {   #wrapper fcn for g(A|W)	
	  	assign(".f_gAW", get(fcn_name))
	  	args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
	  	args <- c(args0, ...)
	  	print(paste("Args passed to", fcn_name, ":"))
	  	print(names(args))
		A <- rbinom(n, 1, do.call(.f_gAW, args))  	
	}		
  	#get all friends Ws as a list of vectors
  	.f.cA_Pa <-function(W1, Net) sapply(Net, function(netwk) W1[netwk], simplify=F)   
  	.f.redefineCov <- function(Var, VarNm) {
		  	#get all friends Ws as a matrix of dim(n,k) filling unused cols with zeros
		  	.f.netCovar <-function(Covar, Net) sapply(Net, function(netwk) c(Covar[netwk],rep(0,k-length(netwk)))) 
			netVar_full <- NULL
			netVar_names <- NULL
			sumnet_name <- NULL
			if (k>0) {
				netVar_full <- .f.netCovar(Var, Net_vect) 
				#(11/05/12) removed the sum of all covariates
				# sumnet_name <- paste("sum_net", VarNm, sep="")
				if (k>1) netVar_full <- t(netVar_full)  #need to transpose matrix when col>1, need a better fix
				netVarNm <- paste("net", VarNm, "_", sep="")
				netVar_names <- paste(netVarNm, c(1:k), sep = "") 
			}
			Var_names <- c(VarNm, netVar_names)	
			d <- cbind(Var, netVar_full)			
			#(11/05/12) removed the sum of all covariates, was messing up h_estimation, same in other fcns
			# Var_names <- c(VarNm, netVar_names, sumnet_name)	
			# if (k>0) d <- cbind(Var, netVar_full, rowSums(netVar_full))
			# else d <- cbind(Var, netVar_full)
			colnames(d) <- Var_names
			return(d)   }
						
  	#Calculate c^Y(Pa(Y_i)) - a function into R^(k+1) that doesn't depend on n or i and is permutation invariant
  	#for each i, get the W's and A's in the network => cY is a vector of (W,A), not including individual (W_i,A_i)
  	.f.cY_Pa <-function(W1, A, Net) {
  			cY_Pa <- list(Pa_Ws=sapply(Net, function(netwk) W1[netwk], simplify=F), Pa_As=sapply(Net, function(netwk) A[netwk], simplify=F)) }  	
  	#Define Y, using the network N, W , A, W1_i & A_i 
	.f.Y <- function(W1, A, cY_Pa) {  	
		  	f.netY_AorWs <- function(netwk, coeff)  sum(coeff * netwk) 		
		  	sum_friendY_Ws <- sapply(cY_Pa$Pa_Ws, f.netY_AorWs, coeff_Y_Wfriend)
		  	sum_friendY_As <- sapply(cY_Pa$Pa_As, f.netY_AorWs, coeff_Y_Afriend)
		  	indivY_W <- coeff_Y_W * W1
		  	indivY_A <- coeff_Y_A * A
		 	Y <- rbinom(n, 1, plogis(Intercept_Y + sum_friendY_Ws + sum_friendY_As + indivY_W + indivY_A))
		 	#print(Intercept_Y + sum_friendY_Ws + sum_friendY_As + indivY_W + indivY_A)
		 	return(Y) 
	}	
  #-----------------------------------------------
  #Generating covars
  #-----------------------------------------------  
  	print("---------------------------")
  	print(paste("f.g_name:", f.g_name))
  	print(paste("f.g_args:", f.g_args))
  	W1 <- .f.W1(EC, n)
   	# print(W1)
  	# Set samp=F to assign same # of friends to each person
  	nFriends <- .f.Net_num(samp=NONFIXED_N_FRIENDS, n, k)   #Generate # of friends to sample for each individ., |N_i|
  	# print(nFriends)
  	ConnectMatx <- .f.genConnectMatx(nFriends)   #Generate the connectivity matrix of friends, N_i
  	# print(ConnectMatx[,])  #print the entire connectivity matrix
  	nFriends <- .f.Net_num_update(ConnectMatx)  #Update # of actual friends sampled, |N_i|
  	# print(nFriends)  	
  	# print(nFriends[which(nFriends!= nFriends_new)] - nFriends_new[which(nFriends!= nFriends_new)])
  	Net_vect <- .f.Net_vect_big(ConnectMatx)    #Get list of vectors of friend's ids for each i
  	# print(Net_vect)  	  	
  	cA_Pa <- .f.cA_Pa(W1, Net_vect)   #Get c^A - fcn for parents of A_i: (W)
  	# print(cA_Pa)  	  	  
  	  
  	if (is.null(f.g_name)) {  
	  	f.g_name <- "f.A"	
	  	A <- .f_g_wrapper(f.g_name) } 
	else { 
	  	A <- .f_g_wrapper(f.g_name, f.g_args)
	}
	#convert f.g_args to text format (for output)
  	if (is.null(f.g_args)) {
	  	f.g_args_txt <- "NA" } 
	else {
	  		f.g_args_txt <- paste(f.g_args, collapse=" ") 
  	} 
  	#Get c^Y fcn for parents of Y_i: (A,W)   
  	cY_Pa <- .f.cY_Pa(W1, A, Net_vect)   
  	W1_netW1 <- .f.redefineCov(W1, "W1")
  	A_netA <- .f.redefineCov(A, "A")

  	Y <- .f.Y(W1, A, cY_Pa)   
 	#Convert N_i, cA_Pa_i, cY_Pa_i to strings
  	Net_str <- .f.mkstrNet(Net_vect)   
  	netW1_str <- .f.mkstrNet(cA_Pa)
  	netA_str <- .f.mkstrNet(cY_Pa$Pa_As)
  	d <- data.frame(C_j, f.g_name, f.g_args_txt, EC, Y, nFriends, W1_netW1, A_netA, Net_str, netW1_str, netA_str, stringsAsFactors = FALSE)
  	return(d) 
  }
  
#Sample the population (several communities) and combine in one dt.frame
gendata_pop <- function(nC = 1, n_arr, k_arr, EC_arr, f.g_list=NULL, f.g_args_list=NULL) {	
	  #n_arr - number of individuals /per community
	  #k_arr - max size of each individual's network /per community
	  #EC_arr - community-level covariate, only influences W_i for i\in nC_j /per community  
	  #nC - # communities to sample
	  .f.EC <- function(i, samp=FALSE, EC_arr) if (samp) sample(EC_arr, 1) else EC_arr[C_j]   
	  EC_rand <- sapply(1:nC, .f.EC, T, EC_arr)  

	  #sample a random assignment of community covars
	  if ((nC>1) & (!(is.null(f.g_list)))) {
		  	fcn_ids <- sample(1:length(f.g_list), nC, replace=TRUE)
		  	f.g_rand <- f.g_list[fcn_ids]
		  	f.g_args_rand <-f.g_args_list[fcn_ids] }
	  else {
	  		if (!(is.null(f.g_list))) {
			  	f.g_rand <- f.g_list
			  	f.g_args_rand <-f.g_args_list }
	  		else {
			  	f.g_rand <- list(NULL)
			  	f.g_args_rand <- list(NULL) }
			}	  
	  if (is.null(f.g_args_list)) f.g_args_rand<-list(NULL)  
	  d <- do.call("rbind", mapply(gendata_Cj, 1:nC, n_arr, k_arr, EC_rand, f.g_rand, f.g_args_rand, SIMPLIFY=FALSE))
	  return(d)  
  	}
  #write.table(d, file = "./netwk_dt.csv", sep = ",")
    

#-------------------------------------------------------------------------------------------------------
# G-Comp & TMLE: Use Monte-Carlo to estimate psi under stochastic intervention g* 
#-------------------------------------------------------------------------------------------------------
  # TO DO LIST:
  #*** 1) Put constraints on h_bar > 0 and < \infty 
  # 3) Check the rate of convergence of Var and Bias (should be 1/n);
  # 4) Allow for W_i's to be independent but not identical (putting mass 1 over each W_i);
#-------------------------------------------------------------------------------------------------------  
  # For given data, take Q[Y|cY]=m.Q.init and calcualte est. of psi under g*=f.g_name using Monte-Carlo integration:
  # * Assume that W_i are iid;
  # * Draw from the distributions of W1 and g*(A|W), keeping N_i constant, recalculate cY and cA each time;
  # * Recalculate Y^c under new W1 and A;
  # * Repeat nrep times and average.
#-------------------------------------------------------------------------------------------------------
  get.MCS_ests <- function(glob.table.h, data, max.err_eps, m.Q.init, m.Q.update_0, m.Q.update, m.gN, nreps.h, n, k, f.g_name, f.g_args, family="binomial") 
  {
		# List of covars (W1 or A) from given network
		.f.get_Pa <-function(Var, Net) sapply(Net, function(netwk) Var[netwk], simplify=F)   
		# Netwk ids strings to lists of vectors
	    .f.mkvecNet <- function(Net_str) lapply(Net_str, function(Net_str_i) as.numeric(unlist(strsplit(Net_str_i, ' ', fixed=TRUE))))
	    
		.f.allCovars <- function(Var, VarNm) { 
			netVar_names <- NULL
			netVar_full <- NULL
			sumnet_name <- NULL
			#mtx <- matrix(0, nrow=n, ncol = k+1) -> try changing to mtx based approach
			if (k>0) {
				netVar_names <- paste(paste("net", VarNm, "_", sep=""), c(1:k), sep = "") 
				# sumnet_name <- paste("sum_net", VarNm, sep="")
				netVar_full <- sapply(Net_vect, function(netwk) c(Var[netwk], rep(0,k-length(netwk))))	
				if (k>1) netVar_full <- t(netVar_full)  #need to transpose matrix when col>1, needs a better fix
			}
			Var_names <- c(VarNm, netVar_names)
			d <- cbind(Var, netVar_full)
			#(11/05/12) remove sum of covars
			# Var_names <- c(VarNm, netVar_names, sumnet_name)	
			# if (k>0) d <- cbind(Var, netVar_full, rowSums(netVar_full))
			# else d <- cbind(Var, netVar_full)
			colnames(d) <- Var_names
			return(d)	
		}
		
	  	# Sample A based on g*=fcn_name, given individ W1's;
	 	.f.gen.A.star <- function(resamp_W1, fcn_name, f_args =NULL) { 
	 		.f_g_wrapper <- function(fcn_name, W1, cA_Pa, n, nFriends, ...) {	
		  		assign(".f_gAW", get(fcn_name))
		  		args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
		  		args <- c(args0, ...)
				do.call(.f_gAW, args) 	}	
	  		resamp_cA_Pa <- .f.get_Pa(resamp_W1, Net_vect)
	  		rbinom(n, 1, .f_g_wrapper(fcn_name, resamp_W1, resamp_cA_Pa, n, data$nFriends, f_args)) 	}
	  		
	 	.f.gen.probA.star <- function(resamp_W1, fcn_name, f_args=NULL) { 
	 		.f_g_wrapper <- function(fcn_name, W1, cA_Pa, n, nFriends, ...) {	
		  		assign(".f_gAW", get(fcn_name))
		  		args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
		  		args <- c(args0, ...)
				do.call(.f_gAW, args)  }	
	  		resamp_cA_Pa <- .f.get_Pa(resamp_W1, Net_vect) 
	  		.f_g_wrapper(fcn_name, resamp_W1, resamp_cA_Pa, n, data$nFriends, f_args)  	}

	  	# Sample A based on estimated regression model for g*;
	 	# .f.gen.A_N <- function(resamp_W1) { 
 			# covars_g <- data.frame(.f.allCovars(resamp_W1, "W1") , nFriends=data$nFriends)	
	  		# return(rbinom(n, 1, predict.glm(m.gN, newdata=covars_g, type="response")))		}
	 	# .f.gen.probA_N <- function(resamp_W1) { 
 			# covars_g <- data.frame(.f.allCovars(resamp_W1, "W1") , nFriends=data$nFriends)			  		
  			# # predA_man <- (expit(model$coefficients %*% rbind(rep(1,n), t(covars_g))))
  			# predict.glm(m.gN, newdata=covars_g, type="response")	}
  			
  		# True prob of Y, given the network, W , A
		.f.get.trueQ0 <- function(samp_data) {  	
				cumsum.matrix=function(datavars, coeff) {
					y <- matrix(1,nrow=dim(datavars)[1],ncol=dim(datavars)[2])
					y[,1] <- datavars[,1]*coeff
					if (dim(datavars)[2] > 1) {
						for (i in 2:dim(datavars)[2]) {
							y[,i] <- y[,i-1] + datavars[,i]*coeff
						}
					}
					return(y[,dim(datavars)[2]])	
				}			
				netWs <- subset(samp_data, select = netW1_1:eval(parse(text=paste("netW1_", k, sep = "")))) 
			  	sum_friendY_Ws <- cumsum.matrix(netWs, coeff_Y_Wfriend)
				netAs <- subset(samp_data, select = netA_1:eval(parse(text=paste("netA_", k, sep = "")))) 
			  	sum_friendY_As <- cumsum.matrix(netAs, coeff_Y_Afriend)
			  	indivY_W <- coeff_Y_W * samp_data$W1
			  	indivY_A <- coeff_Y_A * samp_data$A  		
			 	Y <- plogis(Intercept_Y + sum_friendY_Ws + sum_friendY_As + indivY_W + indivY_A)
			 	return(Y) 
		}

		# Predict probY=Q[Y|A,W] based on g*, W1 and net, for [rep_col] vector of W1
	  	.f.gen.Y <- function(rep_col) { 	
			# Get all friend's W's or A's as a matrix of dim(n,k) filling unused cols with zeros		
	  		samp_dataW1 <- .f.allCovars(resamp_W1_mtx[,rep_col],"W1")  # replace old Ws with resampled values 		  
			samp_dataA <- .f.allCovars(resamp_A_mtx[,rep_col], "A")
			samp_data <- data.frame(C_j=data$C_j, EC=data$EC, nFriends=data$nFriends, samp_dataW1, samp_dataA)
			#-------------------------------------------
			# True Q_0 
	  			TrueQY0 <- .f.get.trueQ0(samp_data)
			#-------------------------------------------
			# G-Comp estimator (est probY based on model for Q_Y)
	  			QY.init <- predict(m.Q.init, newdata=samp_data, type="response")  
			#-------------------------------------------
			# IPTW estimator (est Y_g_star based on model for g_star(A|W)/g(A|W) )
			#........NOT IMPLEMENTED YET.................
			#-------------------------------------------			
			# iid TMLE estimator (under k=0 with h_bar=g*_0/g_0) - 
			# !!!ONLY FOR VERIFICATION WITH NON-MC iid TMLE
			#****???? WHY DOES offset= QY.init on non-logit scale works for iid but not for networkTMLE?
				off <- QY.init
				probA.star <- .f.gen.probA.star(resamp_W1_mtx[,rep_col], f.g_name, f.g_args)
				g.star_0 <- samp_data$A * probA.star + (1-samp_data$A) * (1-probA.star)
				probA <- .f.gen.probA.star(resamp_W1_mtx[,rep_col], "f.A")
				g_0 <- samp_data$A * probA + (1-samp_data$A) * (1-probA)		
				h0 <- g.star_0 / g_0
				if (!is.na(coef(m.Q.update_0))) TMLE_iid <- expit(off + coef(m.Q.update_0) * h0)
			#-------------------------------------------			
			# Network TMLE estimator (adjusted by h_bar ratio and epsilon)	
				if (k>0) {
					cY_samp.mtx <- cbind(ID=c(1: n), 
						subset(samp_data, select = W1:eval(parse(text=paste("netA_", k, sep = ""))))) 
						cY_samp.mtx <- cbind(cY_samp.mtx, nFriends=data$nFriends)
					}
				else {
					cY_samp.mtx <- cbind(ID=c(1: n), subset(samp_data, select = W1:A)) 
					} 
				# replace w/ intersect with local copy  of glob.table.h
				h_bar_list <- est.hbars(glob.table.h, data, cY_samp.mtx, m.gN, nreps.h, n, k, f.g_name, f.g_args)
				h_bars <- h_bar_list$h_bars
				glob.table.h <<- h_bar_list$glob.table.h				
				h <- h_bars$h
				# off_old <- QY.init  #12/19/12 - WRONG, need to covert offset to logistic skale.

				#12/19/12 : 
				# qlogis(p) is the same as the well known ‘logit’ function, 
				# logit(p) = log(p/(1-p)), and plogis(x) has consequently been called the ‘inverse logit’.
				off <- qlogis(QY.init)

				# exp(x)/(1+\exp(x)):
				# if (!is.na(coef(m.Q.update))) QY.star_old <- expit(off + coef(m.Q.update)*h) 
				# plogis: F(x) = 1 / (1 + exp(-(x-m)/s))
				if (!is.na(coef(m.Q.update))) TMLE_net  <- plogis(off + coef(m.Q.update)*h)				
		  	#-------------------------------------------
			return(cbind(MCS.TrueQY0=mean(TrueQY0), MCS.G_comp_net=mean(QY.init), 
						MCS.TMLE_net=mean(TMLE_net), MCS.TMLE_iid =mean(TMLE_iid) )
						)
					# h.star_c=mean(h_bars$h.star_c), h_c=mean(h_bars$h_c), h=mean(h_bars$h)))
					
	  	} # end of .f.gen.Y()
	  	
		#--------------------------------------------------------------------------------------------------	  	
		# Main body of a fcn get.MCS_ests(): MC evalution of the estimators
		#--------------------------------------------------------------------------------------------------	  	
		data <- data[, (names(data) %in% c('C_j','EC','W1','A','nFriends','Net_str'))]
		Net_vect <- .f.mkvecNet(data$Net_str) # get lists of network ids from network strings for each i	
		# Allow this part to loop, until required MCS prob_epsilon for all estimators is reached:
		nrepeat <- 1
		nreps <- 100
		psis_est <- NULL
		repeat {
			# n*nreps matrix of resampled W1s w/ replacement
			resamp_W1_mtx <- replicate(nreps, sample(data$W1, n, replace=TRUE))
			# Pass W1s, one column at a time, resample A's w/ g*
			resamp_A_mtx <- apply(resamp_W1_mtx, 2, .f.gen.A.star, f.g_name, f.g_args)
			
			psis_est <- rbind(psis_est, t(vapply(seq(nreps), .f.gen.Y, 
											c(MCS.TrueQY0=0, MCS.G_comp_net=0, 
											MCS.TMLE_net =0, MCS.TMLE_iid =0))))
										#	,h.star_c=0, h_c=0, h_ratio=0))))
										
			psi_est <- apply(psis_est, 2, mean, na.rm = T)
			var_est <- apply(psis_est, 2, var, na.rm = T) 
			prob_epsilon <- var_est / ((nreps*nrepeat) * (max.err_eps)^2)
			# print(psi_est)			
			# print(prob_epsilon)			
			nrepeat <- nrepeat + 1
			if ( (all(prob_epsilon[c(1:4)] < 0.05)) | (nrepeat >= 50)) {
				print("break reached")
				print(nrepeat)					
				break
			}
		}		
		return(psi_est)
	}
	  	  	
#-------------------------------------------------------------------------------------------------------	
# Estimate h_bar for g_0 (need to be replaced with g_N) and g* given observed data and vector of c^Y's
#-------------------------------------------------------------------------------------------------------	
# METHOD 1 (empirical distribution of \bar{h}=\sum{h_i}):
  # For given data, estimate g[A|cA] and calculate est. of h_i(c) for given value of c and each i. 
  # * Draw B samples from emp distribution of cY, for each i;
  # * Assume the same dimensionality for cY acorss all i;
  # * Assume W1 are iid, use g_N(A|W) and keep N_i constant;
  # * Drawing from distribution of g_N or g* to get A's;
  # * Calculate h_bar from the joint distribution of all c_i^Y over all i;
# METHOD 2:
  # * DESCRIBE ALGORITHM HERE;
#-------------------------------------------------------------------------------------------------------  
  est.hbars <- function(glob.table.h, data, cY.mtx, m.gN, nreps.h, n, k, f.g_name, f.g_args, family="binomial") {
		# List of covars (W1 or A) from given network
	   	.f.get_Pa <-function(Var, Net) sapply(Net, function(netwk) Var[netwk], simplify=F)  

	    # Convert netwk id strings to lists of vectors (of ids)
	    .f.mkvecNet <- function(Net_str) 
	    	lapply(Net_str, function(Net_str_i) as.numeric(unlist(strsplit(Net_str_i, ' ', fixed=TRUE))))

		# Get the network and put all covars in one multidim array (A's or W's), 
		# dim=(nreps , n, k+1)
		.f.arrayCovars <- function(arrVar, Net_vect) { 
				arrVarNet <- array(data=0, dim=c(nreps.h , n, k+1))
				arrVarNet[, , 1] <- arrVar
				for (i in (1:n)) {
					netwk <- Net_vect[[i]]
					if (length(netwk) > 0) {
						arrVarNet[, i, c(2:(length(netwk)+1))] <- arrVar[, netwk]						
					}
				}
				return(arrVarNet)	
		}
		# Get the network and put all covars in one matrix (A's or W's)
		.f.allCovars <- function(Var, VarNm, Net_vect) {
				netVar_names <- NULL
				netVar_full <- NULL
				sumnet_name <- NULL
				#mtx <- matrix(0, nrow=n, ncol = k+1) -> try changing to mtx based approach
				if (k>0) {
					netVar_names <- paste(paste("net", VarNm, "_", sep=""), c(1:k), sep = "") 
					# sumnet_name <- paste("sum_net", VarNm, sep="")
					netVar_full <- sapply(Net_vect, function(netwk) c(Var[netwk], rep(0,k-length(netwk))))	
					if (k>1) netVar_full <- t(netVar_full)  #need to transpose matrix when col>1, needs a better fix
				}
				Var_names <- c(VarNm, netVar_names)
				d <- cbind(Var, netVar_full)
				# (11/05/12) remove sum of covars
				# Var_names <- c(VarNm, netVar_names, sumnet_name)	
				# if (k>0) d <- cbind(Var, netVar_full, rowSums(netVar_full))
				# else d <- cbind(Var, netVar_full)
				colnames(d) <- Var_names
				return(d)	
		}

	  	# Sample A based on g*, given individ W1's
	 	.f.gen.A.star <- function(resamp_W1, Net_vect, fcn_name, f_args =NULL) { 
	 		.f_g_wrapper <- function(fcn_name, W1, cA_Pa, n, nFriends, ...) {	
		  		assign(".f_gAW", get(fcn_name))
		  		args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
		  		args <- c(args0, ...)
				do.call(.f_gAW, args) 	}	
	  		resamp_cA_Pa <- .f.get_Pa(resamp_W1, Net_vect) 
	  		rbinom(n, 1, .f_g_wrapper(fcn_name, resamp_W1, resamp_cA_Pa, n, data$nFriends, f_args)) 	
	  	}
	  		
	 	.f.gen.probA.star <- function(resamp_W1, Net_vect, fcn_name, f_args=NULL) { 
	 		.f_g_wrapper <- function(fcn_name, W1, cA_Pa, n, nFriends, ...) {	
		  		assign(".f_gAW", get(fcn_name))
		  		args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
		  		args <- c(args0, ...)
				do.call(.f_gAW, args)  }	
	  		resamp_cA_Pa <- .f.get_Pa(resamp_W1, Net_vect) 
	  		# print(head(resamp_cA_Pa))
	  		.f_g_wrapper(fcn_name, resamp_W1, resamp_cA_Pa, n, data$nFriends, f_args)  	
	  	}
	  		
	 	# .f.gen.A_N <- function(resamp_W1) { 
 			# covars_g <- data.frame(.f.allCovars(resamp_W1, "W1") , nFriends=data$nFriends)	
	  		# return(rbinom(n, 1, predict.glm(m.gN, newdata=covars_g, type="response")))		}	  		
	 	# .f.gen.probA_N <- function(resamp_W1) { 
 			# covars_g <- data.frame(.f.allCovars(resamp_W1, "W1") , nFriends=data$nFriends)			  		
  			# # predA_man <- (expit(model$coefficients %*% rbind(rep(1,n), t(covars_g))))
  			# predict.glm(m.gN, newdata=covars_g, type="response")	}

		#--------------------------------------------------------------------------------------------------	  		
  		# Matrix form for all g*,g_0=f.A() functions
  		# dim(arr_W1) = (nreps, n, k+1) - array of c_i^Y's over all replications and i=1,..,n
	 	.f.gen.arr.probA.star <- function(arr_W1, Net_vect, NetIdsfull, fcn_name, f_args=NULL) { 
			# Deterministically setting every A=0
			f.A_0 <- function(arr_W1, NetIdsfull, ...) 
						array(data=0, dim=(c(dim(arr_W1)[1], length(NetIdsfull))))	
			# Deterministically setting every A=1
			f.A_1 <- function(arr_W1, NetIdsfull, ...) 
						array(data=1, dim=(c(dim(arr_W1)[1], length(NetIdsfull))))	
			# Usual f.A, g_0		
		 	f.A <- function(arr_W1, NetIdsfull, ...) {		
				k <- dim(arr_W1)[3]-1
				# Get W1's for everyone in F_i (including i)
				sumcoeffW <- coeff_A_W * arr_W1[ , NetIdsfull , 1]
				for (i_k in (1:k)) {
					# Get netW1's for everyone in F_i (including i)
					sumcoeffW <- sumcoeffW + arr_W1[ , NetIdsfull , (i_k+1)]*coeff_A_Wfriend
				}
				x <- Intercept_A + sumcoeffW
		  		pA.mat <- 1 / (1+exp(-(x)))		  				 		
				return(pA.mat)
			}	
			# Stocastic g*, set x% to A=1
			f.A_x <- function(arr_W1, NetIdsfull, x, ...) {
				x <- array(data=x, dim=(c(dim(arr_W1)[1], length(NetIdsfull))))
				return(x)
			}
	 		.f_g_wrapper <- function(arr_W1, NetIdsfull, n, k, nFriends, fcn_name, ...) {	
		  		assign(".f_gAW", get(fcn_name))
		  		args0 <- list(arr_W1=arr_W1, NetIdsfull=NetIdsfull, n=n, k=k, nFriends=nFriends) 
		  		args <- c(args0, ...)
				do.call(.f_gAW, args)  
			}	
	  		.f_g_wrapper(arr_W1, NetIdsfull, n, k, data$nFriends, fcn_name, f_args)  	
	  	}		

		#--------------------------------------------------------------------------------------------------------- 
	  	# Calculate h_bar (under g or g*)	  			
	  	#--------------------------------------------------------------------------------------------------------- 
		#------------------------------------------------------------------------------------------------------				
		# METHOD 2: (11/05/12) New method for estimation of h based on prob of A, see Mark's mtgs notes
		# 12/19/12: VERIFIED THIS WORKS
	  	.f.h.probA <- function(cY.mtx, f.g_name, f.g_args=NULL) {
	  		# Calculate the joint probability of observing a=(a_1,..,a_k) from c^Y_i
	  		# pass a matrix of probA (nreps x k) and fixed vector of (a_1,..,a_k)
			cumprod.matrix <- function(data.indA, data.probA) {
				y <- matrix(1, nrow=dim(data.probA)[1], ncol=dim(data.probA)[2])
				#REPLACE WITH BINOMIAL DENSITY:
				# y[, 1] <- data.probA[,1]*as.integer(data.indA[1]) + 
							# (1-data.probA[,1])*(1-as.integer(data.indA[1]))			
				y[, 1] <- data.probA[,1]^as.integer(data.indA[1]) *
							(1-data.probA[,1])^(1-as.integer(data.indA[1]))										
				if (dim(data.probA)[2] > 1) {
					for (i in 2:dim(data.probA)[2]) {
						#REPLACE WITH BINOMIAL DENSITY:
						# y[,i] <- y[,i-1]*(data.probA[,i]*as.integer(data.indA[i]) + 
											# (1-data.probA[,i])*(1-as.integer(data.indA[i])))	
						y[,i] <- y[,i-1] * (data.probA[,i]^as.integer(data.indA[i]) * 
											(1-data.probA[,i])^(1-as.integer(data.indA[i])))
					}
				}
				return(round(y[,dim(data.probA)[2]],6))	
			}							
			# Get the vectors of network ids
			Net_vect <- .f.mkvecNet(data$Net_str) 	
			# Est probability of W_1=1			
			pW1_1 <- sum(data$W1)/n

			# Split cY.mtx into W and A components
			cY_vec.subsW <- subset(cY.mtx, select = W1:eval(parse(text=paste("netW1_", k, sep = ""))))	
			cY_vec.subsA <- subset(cY.mtx, select = A:eval(parse(text=paste("netA_", k, sep = ""))))	  					
			
			# initialize the vector of h_c values
			h.cY.vec <- NULL
			# loop over cY_vec (c), calculate the prob of observing cY_vec.subsA, 
			# while fixing W's in cY_vec.subsW and permuting all other W's, then sum over i=1,..,N 	
			
			#evaluate h=sum(h_i) for each given c
			.feval_h_c <- function(indA.i) 	{
				indW <- cY_vec.subsW[indA.i,]					
				# use true network W1's:
				c_nFriends <- cY.mtx[indA.i,]$nFriends + 1
				indW_cut <- indW[1:c_nFriends]
				# likelihood of observed W1's								
				p_indW <- (pW1_1)^(sum(indW_cut))*(1-pW1_1)^(sum(1-indW_cut))				
				indA <- cY_vec.subsA[indA.i,]
				# Repeat over all i's, 1,..,N, for which |F_i|= cY.mtx$nFriends			
					# **12/18/12: To speed up try selecting ALL (W1, netW1) for which nFriends=c$nFriends 
					# then creating one array and using it to calculate all h_i's at once for each c
					# -> will eliminate sapply loop below
				.feval_hi <- function(i) 	{
						h_i <- 0
						if (cY.mtx[indA.i,]$nFriends==data[i,]$nFriends) 	{
							NetIdsfull <- c(i, Net_vect[[i]])
							# Assign W component of ciY (indW) to network of i (including i himself)
							samp_W1_mtx.tmp	<- samp_W1_mtx			
							# Replace W1s with c only for j\inF_i (all friends of i and i himself)
							samp_W1_mtx.tmp[, NetIdsfull] <- 
										rep(as.numeric(indW[c(1:length(NetIdsfull))]), each=NROW(samp_W1_mtx))													
							# Construct multidim array (nreps, n, k+1) of all W1's and netW1's
							arr_W1_net <- .f.arrayCovars(samp_W1_mtx.tmp, Net_vect)
							# Matrix (nreps x NetIdsfull) of probA for F_i (g(a)), based on samp_W1 and w
							system.time(prob_A_mtx <- 
											as.matrix(
											.f.gen.arr.probA.star(
												arr_W1_net, Net_vect, NetIdsfull, f.g_name, f.g_args)))
							# Evaluate joint likelihood of (a_1,..,a_k) over all rows of samp_W1_mtx
							cumprobA <- cumprod.matrix(indA, prob_A_mtx)
							# Average over all rows to get h_i			
							h_i <- mean(cumprobA)
						}
						return(h_i)
					}					
				# Evaluate each h_i
				h_i_vec <- sapply(c(1:n), .feval_hi)
				# Sum to get \bar{h_c} from individ h_i's
				h_c <- p_indW * sum(h_i_vec)
				return(h_c)			
			}
			
			# Parallel version, loop over c values
		   	h.cY.vec <- mclapply(c(1:nrow(cY_vec.subsA)), .feval_h_c, mc.cores=processors)
		   	h.cY.vec <- unlist(h.cY.vec)
		   	# Non-parallel version:
		   	# h.cY.vec <- sapply(c(1:nrow(cY_vec.subsA)), .feval_h_c)
			# print("h.cY.vec: ")
			# print(h.cY.vec)
			names(h.cY.vec) <- cY.mtx$ID	
			# Bound \bar{h} here:
			# h.cY.vec[h.cY.vec < 0.01] <- 0.01
			return(h.cY.vec)									
		}		
		# End of .f.h.probA()
		#------------------------------------------------------------------------------------------------------	
		
		#------------------------------------------------------------------------------------------------------	
		# Main body of fcn est.hbars()  		
		#------------------------------------------------------------------------------------------------------	
		.f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))	
		data <- data[, (names(data) %in% c('C_j','EC','W1','A','nFriends','Net_str'))]
						
		#****************
		# Create empty data frames instead:
		# data.frame(a = character(0), b = double(0))
		h_bars <- NULL
		h_bars.prior <-NULL
		
		#****************
		# Define global hash table of cYs and delete cY.mtx members that are already in it
		cY.ID <- .f.mkstrNet(cY.mtx[, !(names(cY.mtx) %in% c('ID'))])
		cY.mtx <- data.frame(cY.ID, cY.mtx)
		#****************
		# Intersect cY with already estimated h_bar vals
		cY.mtx_sel <- cY.mtx[!(cY.mtx$cY.ID %in% glob.table.h$cY.ID), ]

		# run h_bar function only if there are still some values of c_i for which we need to estimate h		
		if (nrow(cY.mtx_sel)>0) {
				# print("cY.mtx_sel")
				# print(cY.mtx_sel)
			#*********sample big matrix of all W's (reps by N)
			samp_W1_mtx <- t(replicate(nreps.h, sample(data$W1, n, replace=TRUE)))			
					# Old method of calculating hbar:
					# time1 <- system.time(P.hbar.c <- .f.h.emp(cY.mtx_sel, nreps.h, "f.A", NULL))
					# time2 <- system.time(P.hbar.str.c <- .f.h.emp(cY.mtx_sel, nreps.h, f.g_name, f.g_args))		
				
			# h_bar under g_0
			time1 <- system.time(P.hbar.c <- .f.h.probA(cY.mtx_sel, "f.A", NULL))
			# h_bar under g_star
			time2 <- system.time(P.hbar.str.c <- .f.h.probA(cY.mtx_sel, f.g_name, f.g_args))
			# P.hbar.str.c <- 0
			h_bars <- data.frame(cY.ID=cY.mtx_sel$cY.ID, h.star_c=P.hbar.str.c, h_c=P.hbar.c, 
								h=(P.hbar.str.c / P.hbar.c))	
			print("time to run hbar.c:")
			print(time1)
			print("time to run hbar.star.c:")
			print(time2)
		}	
		
		if (!is.null(glob.table.h)) {
			h_bars.prior <- merge(cY.mtx[, c("ID", "cY.ID")], glob.table.h, by="cY.ID")
			h_bars.prior <- h_bars.prior[order(h_bars.prior$ID),]
			h_bars.prior <- unique(h_bars.prior)
		}
		
		h_bars_ID <- data.frame(ID=cY.mtx_sel$ID, h_bars)
		h_bars.fin <- rbind(h_bars_ID, h_bars.prior)		
		h_bars.fin <- h_bars.fin[order(h_bars.fin$ID),]
				
		glob.table.h_tmp <- rbind(glob.table.h, h_bars)
		glob.table.h <- glob.table.h_tmp
		return(list(h_bars=h_bars.fin, glob.table.h=glob.table.h))
	}	# End of est.hbars()
#------------------------------------------------------------------------------------------------------	

#-------------------------------------------------------------------------------------------------------	
#Estimate h_bar for g_N and g* given observed data and vector of c^Y's, empirical distribution of h
#-------------------------------------------------------------------------------------------------------	
  # For given data, estimate g[A|cA] and calculate est. of h_i(c) for given value of c and each i. 
  # * Draw B samples from emp distribution of cY, for each i;
  # * Assume the same dimensionality for cY acorss all i;
  # * Assume W1 are iid, use g_N(A|W) and keep N_i constant;
  # * Calculate h_bar from the joint distribution of all c_i^Y over all i;
#-------------------------------------------------------------------------------------------------------  
  est.hbars.emp <- function(data, cY.mtx, m.gN, nreps.h, n, k, f.g_name, f.g_args, family="binomial") {
  	
		.f.h.gA <- function(cY.mtx, f.g_name, f.g_args=NULL) {			
			# List of covars (W1 or A) from given network
		   	.f.get_Pa <-function(Var, Net) sapply(Net, function(netwk) Var[netwk])  		   				   	
			# Get the network and put all covars in one matrix (A's or W's)
			.f.allCovars <- function(Var, VarNm, Net_vect) {
					netVar_names <- NULL
					netVar_full <- NULL
					sumnet_name <- NULL
					#mtx <- matrix(0, nrow=n, ncol = k+1) -> try changing to mtx based approach
					if (k>0) {
						netVar_names <- paste(paste("net", VarNm, "_", sep=""), c(1:k), sep = "") 
						# sumnet_name <- paste("sum_net", VarNm, sep="")
						netVar_full <- sapply(Net_vect, function(netwk) c(Var[netwk], rep(0,k-length(netwk))))	
						if (k>1) netVar_full <- t(netVar_full)  #need to transpose matrix when col>1, needs a better fix
					}
					Var_names <- c(VarNm, netVar_names)
					d <- cbind(Var, netVar_full)
					# (11/05/12) remove sum of covars
					# Var_names <- c(VarNm, netVar_names, sumnet_name)	
					# if (k>0) d <- cbind(Var, netVar_full, rowSums(netVar_full))
					# else d <- cbind(Var, netVar_full)
					colnames(d) <- Var_names
					return(d)	
			}
		  	# Sample A based on g*, given individ W1's
		 	.f.gen.A.star <- function(resamp_W1, Net_vect, fcn_name, f_args=NULL) { 
		 		.f_g_wrapper <- function(fcn_name, W1, cA_Pa, n, nFriends, ...) {	
			  		assign(".f_gAW", get(fcn_name))
			  		args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
			  		args <- c(args0, ...)
					do.call(.f_gAW, args) 	}	
		  		resamp_cA_Pa <- .f.get_Pa(resamp_W1, Net_vect) 
		  		rbinom(n, 1, .f_g_wrapper(fcn_name, resamp_W1, resamp_cA_Pa, n, data$nFriends, f_args)) 	
		  	}		   	
		    # Convert netwk ids strings to lists of vectors
		    .f.mkvecNet <- function(Net_str) lapply(Net_str, function(Net_str_i) 
		    										as.numeric(unlist(strsplit(Net_str_i, ' ', fixed=TRUE))))
		 	# .f.gen.A <- function(resamp_W1) { 
		 		# # Get netWs (as vector of lists)
		 		# samp_dataW1 <- .f.allCovars(resamp_W1,"W1")  # get sampled Ws
	 			# samp_data_g <- data.frame(C_j=data$C_j, EC=data$EC, nFriends=data$nFriends, samp_dataW1)
	 			# # print(samp_data_g)  
	 		  	# probA <- predict.glm(m.gN, newdata=samp_data_g, type="response")
	 			# # print(probA)
		  		# estAs <- rbinom(n, 1, probA)
		  		# return(estAs)
		  		# # return(probA)
	  		# }	  
			# Emp distribution based on sampled A's, using join
	    	# Draw from cY(A,W) based on g, W1 and net, for [rep_col] vector of W1
		  	.f.sample.g.cY_join <- function(rep_col) {
		  		samp_dataW1 <- .f.allCovars(resamp_W1_mtx[,rep_col],"W1", Net_vect)  # get all Ws (with network)
		  		samp_dataA <- .f.allCovars(resamp_A_mtx[,rep_col],"A", Net_vect) # get all As (with network)
		  		# Combine to get a draw of cY(A,W) over all i = a draw of h_bar 
		  		samp_cY <- data.frame(samp_dataW1, samp_dataA) 
		    	# Calculate total no. of samples ci^Y (i=1,..,N) matching each c. This is an estimate of h_bar		
				matches <- join(cY.mtx, samp_cY, type = "inner", match = "all")
				matches <- table(matches$ID)
				matches_count <- rep(0, n)
				names(matches_count) <- c(1:n)		
				matches_count[names(matches)] <- matches
				return(matches_count)
			}
			# Emp distribution based on sampled A's, using merge
		  	.f.sample.g.cY_merge <- function(rep_col) {
		  		samp_dataW1 <- .f.allCovars(resamp_W1_mtx[,rep_col],"W1", Net_vect)  # get all Ws (with network)
		  		samp_dataA <- .f.allCovars(resamp_A_mtx[,rep_col],"A", Net_vect) # get all As (with network)
		  		# Combine to get a draw of cY(A,W) over all i = a draw of h_bar 
		  		samp_cY <- data.frame(samp_dataW1, samp_dataA, nFriends=data$nFriends) 
			  		# print("samp_cY")
			  		# print(head(samp_cY))
			  		# print("cY.mtx")
			  		# print(head(cY.mtx))
		    	# Calculate total no. of samples ci^Y (i=1,..,N) matching each c. This is an estimate of h_bar			
				matches <- merge(cY.mtx, samp_cY, all=F)
				matches <- table(matches$ID)			
					# print("matches")
					# print(matches)	
				matches_count <- rep(0, nrow(cY.mtx))
				names(matches_count) <- cY.mtx$ID			
				matches_count[names(matches)] <- matches
					# print("matches_count")
					# print(matches_count)
				return(matches_count)	
				}	
	
			#------------------------------------------------------------------------------------------------------		  		  	
			# Get lists of network ids from network strings for each i
			Net_vect <- .f.mkvecNet(data$Net_str) 
			#n*nreps matrix of resampled W1's (w/ replacement)
			h.cY.vec <- rep(0, nrow(cY.mtx))
			names(h.cY.vec) <- cY.mtx$ID	
				# print("cY.mtx$ID")
				# print(cY.mtx$ID)
				# print("h.cY.vec")
				# print(h.cY.vec)
			nrepeat <- 0
			#repeat until h for all c's are > 0 (have been observed)
			repeat {
				resamp_W1_mtx <- replicate(nreps.h, sample(data$W1, n, replace=TRUE)) 
				resamp_A_mtx <- apply(resamp_W1_mtx, 2, .f.gen.A.star, Net_vect, f.g_name, f.g_args) 
		  		n_matches_reps <- sapply(seq(nreps.h), .f.sample.g.cY_merge)
		  		h.cY.vec <- h.cY.vec + rowSums(n_matches_reps)
					# print("--------")	  		
					# print(h.cY.vec[h.cY.vec < 0.001])
					# print(all(h.cY.vec >= 0.001))
				nrepeat <- nrepeat + 1
				if ( (all(h.cY.vec >= 0.001)) | (nrepeat >= 10)) {
					print("break reached")
					print(nrepeat)					
					# h.cY.vec <- round(h.cY.vec / (nrepeat * nreps), 4)
					h.cY.vec <- h.cY.vec / (nrepeat * nreps.h)
					print("h.cY.vec")
					print(h.cY.vec)
					h.cY.vec[h.cY.vec < 0.001] <- 0.001
					break
				}
			}
			return(h.cY.vec)	
		}	
		
		#------------------------------------------------------------------------------------------------------		  		
		.f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))			
		data <- data[, (names(data) %in% c('C_j','EC','W1','A','nFriends','Net_str'))]
		cY.ID <- .f.mkstrNet(cY.mtx[, !(names(cY.mtx) %in% c('ID'))])
		cY.mtx <- data.frame(cY.ID, cY.mtx)
		# h_bar under actual g 
		# nreps.h - number of iterations to evaluate h under g_N
			time1 <- system.time(P.hbar.c <- .f.h.gA(cY.mtx, "f.A", NULL))
			# h_bar under g_star
			time2 <- system.time(P.hbar.str.c <- .f.h.gA(cY.mtx, f.g_name, f.g_args))
			# P.hbar.str.c <- 0
			h_bars <- data.frame(cY.ID= cY.mtx$cY.ID, h.star_c=P.hbar.str.c, h_c=P.hbar.c, 
								h=(P.hbar.str.c / P.hbar.c))	
		print("time to run hbar.c:")
		print(time1)
		print("time to run hbar.star.c:")
		print(time2)
		return(h_bars)
	}


# Sample all unique c's, evaluate h and evaluate Q_0(c) under h.star = psi_0
calc_allh <- function(n.h_c.iter, d, n, k, m.Q.init, m.gN, f.g_name, f.g_args) {	
		glob.table.h <- NULL
		
		# Sample c's 
		J <- 5000
		nFriends <- sample(c(0:2),J,replace=TRUE)
	  	.fsampleNet <- function(VarNm) {
				netVar_full <- sapply(nFriends, function(nFriend) {
			  				VarVal<-c(sample(c(0:1),nFriend+1,replace=TRUE), rep(0,k-nFriend))
							return(VarVal) }
			  			)			
				netVarNm <- paste("net", VarNm, "_", sep="")
				netVar_names <- paste(netVarNm, c(1:k), sep = "") 
				Var_names <- c(VarNm, netVar_names)	
				d <- t(netVar_full)
				colnames(d) <- Var_names
				return(d)  
		}
	  	# True prob of Y, given c
		.f.get.trueQ0 <- function(samp_data) {  	
				cumsum.matrix=function(datavars, coeff) {
					y <- matrix(1,nrow=dim(datavars)[1],ncol=dim(datavars)[2])
					y[,1] <- datavars[,1]*coeff
					if (dim(datavars)[2] > 1) {
						for (i in 2:dim(datavars)[2]) {
							y[,i] <- y[,i-1] + datavars[,i]*coeff
						}
					}
					return(y[,dim(datavars)[2]])	
				}			
				netWs <- subset(samp_data, select = netW1_1:eval(parse(text=paste("netW1_", k, sep = "")))) 
			  	sum_friendY_Ws <- cumsum.matrix(netWs, coeff_Y_Wfriend)		  	
				netAs <- subset(samp_data, select = netA_1:eval(parse(text=paste("netA_", k, sep = "")))) 
			  	sum_friendY_As <- cumsum.matrix(netAs, coeff_Y_Afriend)		  	
			  	indivY_W <- coeff_Y_W * samp_data$W1
			  	indivY_A <- coeff_Y_A * samp_data$A  		
			 	Y <- plogis(Intercept_Y + sum_friendY_Ws + sum_friendY_As + indivY_W + indivY_A)
			 	return(Y) 
		}			
		c_vec <- cbind(as.data.frame(.fsampleNet("W1")), as.data.frame(.fsampleNet("A")), nFriends=nFriends)
		c_vec_uniq <- unique(c_vec)
		n_c <- nrow(c_vec_uniq)
		print("# of unique c's: ")
		print(n_c)
			
		# Get all possible c's
		cY.mtx <- cbind(ID=c(1: n_c), c_vec_uniq)
			# print("cY.mtx")
			# print(cY.mtx)
		Q_0 <- .f.get.trueQ0(cY.mtx)	
		QY.init <- predict(m.Q.init, newdata=cY.mtx, type="response")  

		# Evalute h_bar by sampling from its emp. distribution
		h.iter.emp <- 1000
		h_bars.emp <- est.hbars.emp(d, cY.mtx, m.gN, h.iter.emp, n, k, f.g_name, f.g_args) 
		glob.table.h.emp <- h_bars.emp
		c_hbar.emp <- cbind(cY.mtx, h_bars.emp, Q_0=Q_0, prodh_Q_0=Q_0*h_bars.emp$h.star_c)
		print("c_hbar.emp: ")	
		print(c_hbar.emp)
		print("Sum over all h.star_c.emp: ")
		print(sum(c_hbar.emp$h.star_c))
		print("Sum over all h_c.emp: ")
		print(sum(c_hbar.emp$h_c))
		# Using true Q_0 to estimate Q_0.star ~ h.star
		Q_star1.emp  <- sum(c_hbar.emp$prodh_Q_0) / n
		print("Q_star1.emp (direct evaluation of integral):")
		print(Q_star1.emp )
			# Q_star2.emp  <- sum(c_hbar.emp$prodh_Q_0) / sum(c_hbar.emp$h.star_c)
			# print("sum, Q_star2.emp  (weighted mean from distr \bar{h}*): ")
			# print(Q_star2.emp )		
		# Using fitted logistic regression to estimate Q_0.star ~ h.star		
		Q_star.reg1.emp  <- sum(c_hbar.emp$h.star_c*QY.init) / n
		print("Q_star.reg1.emp ")		
		print(Q_star.reg1.emp )	
			# Q_star.reg2.emp  <- sum(c_hbar.emp$h.star_c*QY.init) / sum(c_hbar.emp$h.star_c)
			# print("Q_star.reg2.emp ")		
			# print(Q_star.reg2.emp )	
		# Using true Q_0 to estimate Q_0 ~ h (g_0)
		Q_1.emp <- sum(Q_0 * c_hbar.emp$h_c) / n
		print("Q_1.emp (under h, direct evaluation of integral):")
		print(Q_1.emp)
			# Q_2.emp <- sum(Q_0 * c_hbar.emp$h_c) / sum(c_hbar.emp$h_c)
			# print("Q_2.emp (under h, weighted mean from distr \bar{h}): ")
			# print(Q_2.emp)
		# Using fitted logistic regression to estimate to estimate Q_0 ~ h (g_0)  
		Q_reg1.emp <- sum(c_hbar.emp$h_c*QY.init) / n
		print("Q_reg1.emp")		
		print(Q_reg1.emp)	
			# Q_reg2.emp <- sum(c_hbar.emp$h_c*QY.init) / sum(c_hbar.emp$h_c)
			# print("Q_reg2.emp")		
			# print(Q_reg2.emp)			
		
		# Evalute h_bar by calculating prob of A (g(A|W))
		h_bar_list <- est.hbars(glob.table.h, d, cY.mtx, m.gN, n.h_c.iter, n, k, f.g_name, f.g_args)
		glob.table.h <- h_bar_list$glob.table.h
		h_bars <- h_bar_list$h_bars
		# Evaluate  1/N (\barQ(c)h*(c)) 
		c_hbar <- cbind(cY.mtx, h_bars, Q_0=Q_0, prodh_Q_0=Q_0*h_bars$h.star_c)
		print("c_hbar: ")	
		print(c_hbar)
		print("Sum over all h.star_c: ")
		print(sum(c_hbar$h.star_c))
		print("Sum over all h_c: ")
		print(sum(c_hbar$h_c))
		# Using true Q_0 to estimate Q_0.star ~ h.star
		Q_star1 <- sum(c_hbar$prodh_Q_0) / n
		print("Q_star1 (direct evaluation of integral):")
		print(Q_star1)
		# # Q_star2 <- sum(c_hbar$prodh_Q_0) / sum(c_hbar$h.star_c)
		# # print("sum, Q_star2 (weighted mean from distr \bar{h}*): ")
		# # print(Q_star2)
		# Using fitted logistic regression to estimate Q_0.star ~ h.star		
		Q_star.reg1 <- sum(c_hbar$h.star_c*QY.init) / n
		print("Q_star.reg1")		
		print(Q_star.reg1)	
			# Q_star.reg2 <- sum(c_hbar$h.star_c*QY.init) / sum(c_hbar$h.star_c)
			# print("Q_star.reg2")		
			# print(Q_star.reg2)	
		# Using true Q_0 to estimate Q_0 ~ h (g_0)
		Q_1 <- sum(Q_0 * c_hbar$h_c) / n
		print("Q_1 (under h, direct evaluation of integral):")
		print(Q_1)
			# Q_2 <- sum(Q_0 * c_hbar$h_c) / sum(c_hbar$h_c)
			# print("Q_2 (under h, weighted mean from distr \bar{h}): ")
			# print(Q_2)
		# Using predicted Q to estimate Q_0 ~ h (g_0)		
		Q_reg1 <- sum(c_hbar$h_c*QY.init) / n
		print("Q_reg1")		
		print(Q_reg1)	
			# Q_reg2 <- sum(c_hbar$h_c*QY.init) / sum(c_hbar$h_c)
			# print("Q_reg2")		
			# print(Q_reg2)	

		print("h_bars.emp: ")
		print(h_bars.emp)
		
		print("glob.table.h: ")
		print(glob.table.h)
		
		print(failvar)		
		return(list(glob.table.h = glob.table.h, 
					Q_star1 = Q_star1, Q_star.reg1 = Q_star.reg1, 
					Q_1=Q_1, Q_reg1=Q_reg1))
}
	
#-------------------------------------------------------------------------------------------------------	
# TMLE ESTIMATOR
#-------------------------------------------------------------------------------------------------------	
estimate_tmle <- function(psi_0, Qform, gform, d, nC=1, n, k, EC, f.g_name=NULL, f.g_args=NULL, family="binomial") {
		#----------------------------------------------------------------------------------
		# MONTE-CARLO SIMULATION PARAMETERS
		#----------------------------------------------------------------------------------		
		max.err_est <- 0.01	#maximum absolute error for the MCS estimator
		#iterations for computing baseline h
		n.h_c.iter <- 500
		#iterations for computing h for g-comp
		n.reps.Qh <- 500 #h MC evaluation during Q.star update
		#----------------------------------------------------------------------------------
		print("psi_0: ") 				
		print(psi_0) 				
		# reset global h_bar hash table
		#****************
		# Create empty data frames instead:
		# data.frame(a = character(0), b = double(0))
		glob.table.h <<- NULL
		#****************
	   	.f.est <- function(data, form, family) return(glm(as.formula(form), data=data, family=family))
	   	
		#-------------------------------------------					
		# 1) Specify model for Q(Y|A,W)		
    	m.Q.init <- .f.est(d, Qform, family= family)
    	QY.init <- predict(m.Q.init, newdata=d, type="response")
    	# print(summary(m.Q.init))
		d <- data.frame(d, QY.init)
		#-------------------------------------------							
		# 2) Specify model for g(A,W)			
		m.gN <- .f.est(d, gform, family= family)
    	# print(summary(m.gN))
		#-------------------------------------------					
		# 3) Calculate h, h*, h^tilda for all unique values of c
		allh_ojb <- calc_allh(n.h_c.iter, d, n, k, m.Q.init, m.gN, f.g_name, f.g_args)
    	glob.table.h <- allh_ojb$glob.table.h
		#-------------------------------------------	
		# 4a) Estimate standard TMLE (under iid data) with stochastic intervention g*
	 	.f.gen.probA.star <- function(fcn_name, f.g_args=NULL) { 
	    	.f.mkvecNet <- function(Net_str) 
	    		lapply(Net_str, function(Net_str_i) as.numeric(unlist(strsplit(Net_str_i, ' ', fixed=TRUE))))
	 		.f.get_Pa <-function(Var, Net) sapply(Net, function(netwk) Var[netwk], simplify=F)  
	 		.f_g_wrapper <- function(fcn_name, W1, cA_Pa, n, nFriends, ...) {	
		  		assign(".f_gAW", get(fcn_name))
		  		args0 <- list(W1=W1, cA_Pa=cA_Pa, n=n, nFriends=nFriends) 
		  		args <- c(args0, ...)
				do.call(.f_gAW, args)  }
			Net_vect <- .f.mkvecNet(d$Net_str) 
	  		cA_Pa <- .f.get_Pa(d$W1, Net_vect) 
	  		.f_g_wrapper(fcn_name, d$W1, cA_Pa, n, d$nFriends, f.g_args)  	
	  		}	
		probA.star <- .f.gen.probA.star(f.g_name, f.g_args)
		g.star_0 <- d$A * probA.star + (1-d$A) * (1-probA.star)
		probA <- .f.gen.probA.star("f.A")
		g_0 <- d$A * probA + (1-d$A) * (1-probA)		
		h_bars_0 <- data.frame(h.star_c=g.star_0, h_c=g_0, h=(g.star_0 / g_0))
		
		Y <- d$Y
		off <- QY.init 
		h_0 <- h_bars_0$h
		m.Q.update_0 <- glm(Y ~ -1+h_0, data=data.frame(Y, off, h_0), family=family, offset = off)		
		probA.star <- .f.gen.probA.star(f.g_name, f.g_args)
		probA <- .f.gen.probA.star("f.A")
		QY.init_0 <- predict(m.Q.init, newdata=transform(d,A=1), type="response") * probA.star +
						predict(m.Q.init, newdata=transform(d,A=0), type="response") * (1-probA.star)
		QY.star_iid <- (expit(predict(m.Q.init, newdata=transform(d,A=1), type="response") 
					+ coef(m.Q.update_0) * probA.star/probA) * probA.star) + 
						(expit(predict(m.Q.init, newdata=transform(d,A=0), type="response") 
					+ coef(m.Q.update_0) * (1-probA.star)/(1-probA)) * (1-probA.star))					
		# print(paste("Psi gcomp_iid for stochastic intervention:", mean(QY.init_0)))
		# print(paste("Psi TMLE_iid for stochastic intervention:", mean(QY.star_iid)))				
		#-------------------------------------------					
		# 4b) Estimate network TMLE update step (logistic regression on h with offset as QY.init)
		# Estimate h_bar and h*_bar for each c_i
		if (k>0) cY.mtx <- cbind(ID=c(1:n), 
							subset(d, select = W1:eval(parse(text=paste("netA_", k, sep = "")))))
		else cY.mtx <- cbind(ID=c(1: n), subset(d, select = W1:A))		
		cY.mtx <- cbind(cY.mtx, nFriends=d$nFriends)
		
		#*** Replace w/ intersect from local copy of glob.table.h	
		h_bar_list <- est.hbars(glob.table.h, d, cY.mtx, m.gN, n.h_c.iter, n, k, f.g_name, f.g_args) 
		h_bars <- h_bar_list$h_bars
		d <- data.frame(d, h_bars) 
		Y <- d$Y	
		off_old <- QY.init 
		# print("off_old: ")
		# print(off_old)
		off <- qlogis(QY.init)
		# print("off: ")
		# print(off)
		h <- h_bars$h
		# m.Q.update_old <- glm(Y ~ -1 + h, data=data.frame(Y, off, h), family=family, offset = off)
        m.Q.update <- glm(Y ~ -1 + h + offset(off), data=data.frame(Y, off, h), family=family)		
		
		# if (!is.na(coef(m.Q.update))) QY.update_old <- expit(off + coef(m.Q.update)*h) 
		if (!is.na(coef(m.Q.update))) QY.update  <- plogis(off + coef(m.Q.update)*h)
		print("New TMLE fit for h_ratio (espilon_hat): ")
		print(unlist(coef(m.Q.update)))
		d <- data.frame(d, QY.update)	
		# print("d, h_bars")		
		# print(d)				
		#-------------------------------------------									
		# 5) Run Gcomp, TMLE_net estimator with an updated QY.star, recalculating h_bar each time
		MCS_ests <- get.MCS_ests(glob.table.h, d, max.err_est, m.Q.init, m.Q.update_0, 
						m.Q.update, m.gN, n.reps.Qh, n, k, f.g_name, f.g_args, family)
		print("Monte Carlo Sim ests:")
		print(MCS_ests)		
			#-------------------------------------------	
				# *** Verify Double Robustness
				Y <- d$Y
				# psi_0 <- 0.35995
				# psi_0 <- mean(gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)$Y)
				# print("psi_0: ")
				# print(psi_0)					
				DR_eq_Q0 <- mean(h*Y)
				DR_eq_Qmiss <- mean(h*QY.init)				
				print("DR_equation_est_Q0: ")				
				print(DR_eq_Q0)
				print("DR_equation_est_Qmiss:")				
				print(DR_eq_Qmiss)		
				print("Q0-Qmiss: ")				
				print(DR_eq_Q0-DR_eq_Qmiss)		
			#-------------------------------------------	
			
		all_ests <- c(MCS_ests, 
							DR.Q0=allh_ojb$Q_1, DR.Q0.reg=allh_ojb$Q_reg1, 
							DR.Q_star=allh_ojb$Q_star1, DR.Q_star.reg=allh_ojb$Q_star.reg1, 
							DR_eq_Q0 = DR_eq_Q0, DR_eq_Qmiss = DR_eq_Qmiss,
							TMLE_iid_g0=mean(QY.star_iid) )
							
			# g_comp.iid=mean(QY.init_0), tmle_eps_0=coef(m.Q.update_0), tmle_eps=coef(m.Q.update),
			# h_tilda=mean(d$h), hstar=mean(d$h.star_c), h=mean(d$h_c),
		print("All ests:")
		print(all_ests)		
		return(list(all_ests=all_ests, h_bars=h_bars, h_bars_0=h_bars_0))
	}

#-------------------------------------------------------------------------------------------------------
#Calculate psi_0 and repetitevly simulate tmle psi est. under given stochastic intervention g* (g_fcn) 
#-------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
#Assess the performance of the TMLE and G-comp for different g*, N and K
#----------------------------------------------------------------------------------
  run_sim <- function() {
  	#----------------------------------------------------------------------------------
  	# SIMULTAION PARAMETERS
  	#----------------------------------------------------------------------------------
  	set.seed(20)
  	k_pop_arr <- c(2)	  
	# narr <- c(50, 100, 200)
	# n_sims <- 200
	narr <- c(100)
	n_sims <- 1
	
 	# nrep_psi0 <- 5
	nrep_psi0 <- 1   #no. of times repeat calc. of psi_0
  	n_psi0 <- 1000  #pop size for calc of psi_0
	outfolder <- "SimResults"
	resfname <- paste("./", outfolder, "/PsiSimRes_TexTabs.txt", sep="")
	
	args_list = list(list(x=0.2))
	g_fcns_list_sim <- list(f.A_x="f.A_x")
		# args_list = list(list(x=0.2), list(x=0.4), list(x=0.6),  NULL)
		# g_fcns_list_sim <- list(f.A_x="f.A_x", f.A_x="f.A_x", f.A_x="f.A_x", f.A_nospill ="f.A_nospill")
		
		# args_list = list(NULL)
		# g_fcns_list_sim <- list(f.A_nospill ="f.A_nospill")
		
		# args_list = list(NULL)
		# g_fcns_list_sim <- list(f.A_0="f.A_0")

	f.g_args_0 =list(NULL)
	f.g_name_0 = "f.A"

  	#----------------------------------------------------------------------------------
	# Specify Qform for each i
  	qform_Vars_arr <- c("W1 + A", "A")
   	# qform_Vars_arr <- c("W1 + A")
  	# qform_Vars_arr <- c("A")		
	# Define Qform here
	.f.Qform <- function(k, qform_Vars) {
		qform_NetVars <- NULL
		if (k > 0) {
			#(11/05/12) removed Qform based on sum of covars
			# qform_NetVars <- paste("sum_netW1", "sum_netA", "nFriends", sep = " + ")
			qform_NetVars <- paste(c(paste("netW1_", c(1: k), sep = ""), 
			paste("netA_", c(1: k), sep = ""), "nFriends"), collapse=" + ")	
		}
		qform_full <- paste(c(qform_Vars, qform_NetVars), collapse = " + ")
		Qform <- paste("Y ~ ", qform_full, collapse = "")
		print(Qform)	
		return(Qform)	
	}
	# Define gform here
	.f.gform <- function(k, gform_Vars){
		gform_NetVars <- NULL
		if (k > 0) {
			#(11/05/12) removed gform based on sum of covars
			# gform_NetVars <- paste("sum_netW1", "nFriends", sep = " + ")
			gform_NetVars <- paste(c(paste("netW1_", c(1: k), sep = ""), "nFriends"), collapse=" + ")
		}
		gform_full <- paste(c(gform_Vars, gform_NetVars), collapse = " + ")
		gform <- paste("A ~ ", gform_full, collapse = "")				
		print(gform)
		return(gform)	
	}
  	#----------------------------------------------------------------------------------

  	#----------------------------------------------------------------------------------
  	# DATA SIMULATION
  	#----------------------------------------------------------------------------------		  
	.f.simulate_tmle <- function(n_arr, k_pop, EC_pop, f.g_name, f.g_args=NULL, nC=1, outfolder) {		
	  	# Summary stats for the estimator
		summ <- function(res, psi_0) { 
			.f_calcsumm <- function(ests) {
				psi_N <- mean(ests)
				bias_N <- psi_N - psi_0
				bias_perc <- 100*bias_N/psi_0
				var_N <- var(ests)			
				mse_N <- bias_N^2 + var_N				
				# CoefVar <- 100*sd(res)/mean(res)	
				return(cbind(psi_0=round(psi_0, 4), psi_est=round(psi_N, 4), 
							bias=round(bias_N, 4), bias_perc=round(bias_perc, 4), 
							Var=round(var_N, 4), MSE=round(mse_N, 4)))
			}
			all_ests <- res$all_ests		
			tab_new <- NULL
			for (row in c(1:nrow(all_ests))) {
				 tab_tmp <- .f_calcsumm(all_ests[row, ])
				 tab_new <- rbind(tab_new, tab_tmp)
			}
			rownames(tab_new)  <- rownames(all_ests)
			return(tab_new) 	
		}
	    # Estimate psi n_sims number of times (simulating new dataset each time)
		# repeat_sim <- function(n, n_sims) replicate(n_sims, sim_estimator(n,k_pop,EC_pop,f.g_name,f.g_args))
	 	repeat_sim <- function(n, n_sims, psi_0) {
		 	# Calculate psi by generating a sample data from actual data gen. distrib.
			sim_estimator <- function(n, k, EC, f.g_name_0, f.g_args_0=NULL, f.g_name, f.g_args=NULL, nC=1) {
				d <-gendata_pop(nC, n, k, EC, f.g_name_0, f.g_args_0)
				# print(d)
				Qform <- .f.Qform(k, qform_Vars)
				print(Qform)
				gform_Vars <- "W1"
				gform <- .f.gform(k, gform_Vars)
				print(gform)
				tmle_out <- estimate_tmle(psi_0, Qform, gform, d, nC=1, n, k, EC, f.g_name, f.g_args)
				# print(tmle_out)
				return(tmle_out) 
			}	
			# Get all output from tmle fcn as list for each sim

			# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			# 12/23/12 CHANGED TO RUN SIMULATIONS ONLY ON 1 CORE, NO PARALLEL
			# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		   	tmle_sims <- mclapply(seq(n_sims), 
		   		function(x) sim_estimator(n, k_pop, EC_pop, f.g_name_0, f.g_args_0, f.g_name, f.g_args), mc.cores=1)
			# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
		   	# Comb through output over sims, put into array
			all_ests <- sapply(tmle_sims, function(x) x$all_ests)	
			print("sims all_ests: ")			
			print(all_ests)
			return(list(all_ests=all_ests))	
		}
  		#----------------------------------------------------------------------------------		  
	  	# Calculate the true value of the parameter psi_0 under g* from the true data gen. distribution  	
	  	assign(paste("reps_psi_0_", f.g_name, sep=""), 
	  				unlist(mclapply(seq(nrep_psi0), function(x) 
	  						mean(gendata_pop(nC, n_psi0, k_pop, EC_pop, f.g_list=f.g_name, 
	  										list(f.g_args))$Y), mc.cores=processors)))
	  	assign(paste("psi_0_", f.g_name, sep=""), 
	  				mean(get(paste("reps_psi_0_", f.g_name, sep=""))))
	  	print(get(paste("psi_0_", f.g_name, sep="")))
	  	#----------------------------------------------------------------------------------		  

		# Calculate the the estimators (psi_est object) n_sims times, for each N_pop and qform
	  	res_tab_full <- NULL
	  	res_tab_rnames <-NULL
	  	for (n in n_arr) {
	  		for (qform_Vars in qform_Vars_arr) {
		  		print(paste("repeat_sim time, n=", n)) 
		  		# 12/17/12: REWRITE, THIS IS A VERY DUMB WAY TO DO IT!!!!
		  		print(system.time(assign(paste("res_n", n,"_k",k_pop,"_",f.g_name,sep=""), 
		  								repeat_sim(n, n_sims, get(paste("psi_0_", f.g_name, sep=""))))))  
			  	sink(paste("./", outfolder, "/simres_k_", k_pop, "_", f.g_name, ".txt",sep=""), append=T)	
			  	
			  	#save results of this run
		  		titleTxt <- paste("n = ", n, "; k = ", k_pop, "; Q.Y = ", qform_Vars, "; n_sims = ", n_sims, sep="")
		  		titleTeX <- paste("K =", k_pop, "; N =", n, "; Q.Y = ", qform_Vars)
		  		res_tab_rnames <- c(res_tab_rnames, titleTeX)
		  		
				res_tab <- summ(get(paste("res_n", n, "_k",k_pop,"_",f.g_name,sep="")), 
								get(paste("psi_0_", f.g_name, sep="")))
								
				res_tab_full <- rbind(res_tab_full, res_tab)
				cat("\n", titleTxt, "\n")
				print(res_tab) 		
				sink()

				# Also print to usual output stream
				print(titleTxt)
				print(res_tab) 
	  		}
	  	}	 		  	    	
		# rownames(res_tab_full) <- res_tab_rnames
		# print(res_tab_full)
		# TexRes_tab <-xtable(res_tab_full, digits=4)
		# print.xtable(TexRes_tab, type="latex", file="./SimResults/SimRes_TexTabs.txt", caption.placement="top", 
		# include.rownames=TRUE, append=TRUE)		
		return(res_tab_full)  	
	}
	# end of .f.simulate_tmle() - main simulation runner
  	#----------------------------------------------------------------------------------		  
	
	sink(resfname)
	sink()	

  	#----------------------------------------------------------------------------------		  	
	# Calculate the the estimators (.f.simulate_tmle() fcn) for each g_fcn and each k
	for (g_fcn_i in c(1:length(g_fcns_list_sim))) {
		res_tab_fcn <- NULL
	  	for (k in k_pop_arr) {
	  		print("run sim w/ these g*")
	  		print(g_fcns_list_sim[[g_fcn_i]])
	  		print(args_list[[g_fcn_i]])  	
	  		print(system.time(res_tab <- .f.simulate_tmle(narr, k_pop=k, EC_pop=0.35, 
	  							f.g_name=g_fcns_list_sim[[g_fcn_i]], f.g_args= args_list[[g_fcn_i]], 
	  							nC=1, outfolder)))
	  		res_tab_fcn <- rbind(res_tab_fcn, res_tab)  		
		}
		
		TeXres_tab_fcn <-xtable(res_tab_fcn, digits=4, 
							caption=paste("g* name: ", "$", g_fcns_list_sim[[g_fcn_i]], "$", ", Args: ", 
									paste(args_list[[g_fcn_i]], collapse=" "), ", N sims:", n_sims))
		print.xtable(TeXres_tab_fcn, type="latex", file=resfname, caption.placement="top", 
						include.rownames=TRUE, append=TRUE,  
						hline.after=c(-1,seq(0,length(narr)*length(k_pop_arr), by=length(narr))) )	 
	}
  	#----------------------------------------------------------------------------------		  	
	
  }  

run_sim()
   

#---------------------------------------------------------------------------------
# Examples of counfounding and interference effect in k=2
#---------------------------------------------------------------------------------
g_fcns_list <- list(f.A="f.A", f.A_0="f.A_0", f.A_1="f.A_1", f.A_cutt_offW="f.A_cutt_offW", f.A_nospill="f.A_nospill", f.A_spillonly="f.A_spillonly", f.A_x="f.A_x", f.A_xlevelW="f.A_xlevelW" , f.A_xlevelNi="f.A_xlevelNi", f.A_1_highNi="f.A_1_highNi")	
g_fcns_list_sim <- list(f.A_0="f.A_0", f.A_1="f.A_1", f.A_x="f.A_x", f.A_1_highNi="f.A_1_highNi", f.A_nospill="f.A_nospill")
args_list <- list(f.A_0=NULL, f.A_1 =NULL, f.A_x=list(x=0.2), f.A_1_highNi=NULL, f.A_nospill=NULL)

#---------------------------------------------------------------------------------
# k <- 0
#---------------------------------------------------------------------------------
# mydat <-gendata_pop(1, 20000, k, 0.35)
# mean(mydat$A)
# [1] 0.5484
# mean(mydat$Y)
# [1] 0.47865
# head(mydat)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_x 
# f.g_args <-  args_list$f.A_x
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# mean(gen.gstar$A)
# # [1] 0.201
# mean(gen.gstar$Y)
# # [1] 0.29305
# head(gen.gstar)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_0
# f.g_args <-  args_list$f.A_0
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# > mean(gen.gstar$A)
# [1] 0
# > mean(gen.gstar$Y)
# [1] 0.1835
# head(gen.gstar)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_1 
# f.g_args <-  args_list$f.A_1
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# > mean(gen.gstar$A)
# [1] 1
# > mean(gen.gstar$Y)
# [1] 0.72595
# head(gen.gstar)
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# k <- 2
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A
# f.g_args <- NULL
# mydat <-gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# mean(mydat$A)
# [1] 0.60785
# mean(mydat$Y)
# [1] 0.60535
# [1] 0.60685
# head(mydat)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_0
# f.g_args <-  args_list$f.A_0
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# mean(gen.gstar$A)
# # [1] 0
# mean(gen.gstar$Y)
# # [1] 0.228
# head(gen.gstar)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_1 
# f.g_args <-  args_list$f.A_1
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# mean(gen.gstar$A)
# # [1] 1
# mean(gen.gstar$Y)
# # [1] 0.85995
# head(gen.gstar)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_x 
# f.g_args <-  args_list$f.A_x
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# > mean(gen.gstar$A)
# [1] 0.19555
# > mean(gen.gstar$Y)
# [1] 0.35995
# head(gen.gstar)
#---------------------------------------------------------------------------------
# f.g_name <- g_fcns_list_sim$f.A_nospill 
# f.g_args <-  args_list$f.A_nospill
# gen.gstar <- gendata_pop(1, 20000, k, 0.35, f.g_name, f.g_args)
# > mean(gen.gstar$A)
# [1] 0.5555
# > mean(gen.gstar$Y)
# [1] 0.57015
# head(gen.gstar)

