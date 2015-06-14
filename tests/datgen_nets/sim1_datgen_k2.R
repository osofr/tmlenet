#------------------------------------
# sim_genk2dat.R
# author: Oleg Sofrygin 
# DATA SIMULATOR FOR NETWORKS WITH K=2 (01/02/2014)
# Source the file to generate a network using these functions
#------------------------------------

#----------------------------------------------------------------------------------
#DEFINE VARIOUS REGRESSION FORMULAS FOR Q AND g
#----------------------------------------------------------------------------------
.f.Qform <- function(k, qform_Vars, miss=FALSE) {
	qform_NetVars <- NULL
	if ((k > 0) & (!miss)) {
		qform_NetVars <- paste(c(paste("netW1_", c(1: k), sep = ""),
		paste("netA_", c(1: k), sep = ""), "nFriends"), collapse=" + ")
	}
	return(paste("Y ~ ", paste(c(qform_Vars, qform_NetVars), collapse = " + "),
																		collapse = ""))
}
.f.Qform2 <- function(misscat=1, k, vcovars_i) {
  	# misscat = 1 - correctly specificied model for Q
  	# misscat = 2 - only i's covariates + no. of friends
  	# misscat = 3 - interactions between all treatments and W1's + no. of friends
	v_fullnetcovars <- c(paste("netW1_", c(1: k), sep = ""), paste("netA_", c(1: k), sep = ""))
	# print(v_fullnetcovars)
	i_interact <- paste("I(",vcovars_i[1],"*",vcovars_i[2], ")", sep="")
	# print(i_interact)
	v_fullnetinteract <- paste(rep("I(netW1_",k), c(1:k), "*", rep("netA_",k), c(1: k), ")", sep = "")
	# print(v_fullnetinteract)
	model_Q <- switch(misscat, 
						# 1:
						paste(c(vcovars_i, v_fullnetcovars, "nFriends"),collapse="+"),
						# 2:
						paste(c(vcovars_i, "nFriends"), collapse="+"),
						# 3:
						paste(c(i_interact, v_fullnetinteract, "nFriends"), collapse="+")
					)
	return(paste("Y ~ ", model_Q, sep=""))
}
.f.gform <- function(k, k_miss=1, gform_Vars, miss=F) {
	gform_NetVars <- NULL
	if (miss) k <- k_miss
	if (k > 0) {	
		gform_NetVars <- paste("netW1_", c(1: k), sep = "")
		if (!miss) gform_NetVars <- c(gform_NetVars, "nFriends")
	}
	gform <- paste("A ~ ", paste(c(gform_Vars, gform_NetVars), collapse = " + "),
																collapse = "")
	return(gform)
}
  #----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# DATA GENERATING PARAMETERS
# FUNCTIONS FOR TREATMENT MECHANISM g(A|W) AND g*(A|W)
#---------------------------------------------------------------------------------
# g_0(A_i|W)=P[A=1|cA]:  
f.A <- function(k, data, ...) {
	#---------------------------------------------------------------------------------
  	# coefficients for g(A|W)
  	Intercept_A <- -1.2
  	coeff_A_W <- 1.5
  	# coeff_A_W <- 0.8
  	# for k=2, no PV
  	coeff_A_Wfriend <- 0.6
  	# for k=2, near PV
  	# coeff_A_Wfriend <- 1
  	# for k=6
  	# coeff_A_Wfriend <- 0.20  	
	#---------------------------------------------------------------------------------
  	n <- nrow(data)
  	nFriends <- data$nFriends
  	W1 <- data$W1
	netW1 <- subset(data, select = netW1_1:eval(parse(text=paste("netW1_", k, sep = ""))))
  	.f_evalA_mtx <- function(W1, netW1, n, nFriends, ...) {
  		# print("netW1"); print(head(netW1))
	  	sum_friendWs <- matrix(0, nrow=n, ncol=1)
	  	# sum_friendW2 <- matrix(0, nrow=n, ncol=1)
		for (k_ind in (1:dim(netW1)[2])) {
			sum_friendWs <- sum_friendWs + netW1[, k_ind] * coeff_A_Wfriend				
		}
		return(sum_friendWs)
  	}
  	# P(A) component from individual infection risk (W1)
	indivW <- coeff_A_W * W1
 	sum_friendWs <- .f_evalA_mtx(W1, netW1, n, nFriends, ...)
  	probA <- plogis(Intercept_A + indivW + sum_friendWs)	
 	return(probA)   	
}
# Deterministically set every A=0
f.A_0 <- function(data, ...) rep(0,nrow(data))		
# Deterministically set every A=1
f.A_1 <- function(data, ...) rep(1,nrow(data))
# Set x% of community to A=1 (returns probability P(A=1))
f.A_x <- function(data, x, ...) rep(x, nrow(data)) 
# Set x% of community to A=1 stratified by W_i
f.A_xlevelW1 <- function(data, xW1_0, xW1_1, ...) {
	n <- nrow(data)
	W1 <- data$W1 		
  	x <- rep(0,n)
	x[which(W1==0)] <- xW1_0
	x[which(W1==1)] <- xW1_1
  	return(x)   
}
# Set x% of community to A=1 based on nFriends={0,1,2}
f.A_xlevelnFriends <- function(data, x_nF0, x_nF1, x_nF2, ...) {	
  	n <- nrow(data)
  	nFriends <- data$nFriends	
  	x <- rep(0,n)
	x[which(nFriends == 0)] <- x_nF0
	x[which(nFriends == 1)] <- x_nF1
	x[which(nFriends == 2)] <- x_nF2
  	return(x)
}
# Set x% of community to A=1 based on connectivity |N_i|= {low, high}, based on median
f.A_xlevelNi <- function(data, x_Ni_low, x_Ni_high, ...) {	
  	n <- nrow(data)
  	nFriends <- data$nFriends	
    Ni_med = quantile(nFriends, 0.5)
    #print(paste("Ni_low", Ni_med))
  	x <- rep(0,n)
	x[which(nFriends <= Ni_med)] <- x_Ni_low
	x[which(nFriends > Ni_med)] <- x_Ni_high
  	return(x)  
}
# Deterministically set A=1 based on cutt-off value for Ws or if W_i = 1
f.A_cutt_offW <- function(k, data, cutt_offW, ...) {	
  	n <- nrow(data)
  	W1 <- data$W1
  	netW1 <- subset(data, select = netW1_1:eval(parse(text=paste("netW1_", k, sep = ""))))
  	sum_friendWs <- matrix(0, nrow=n, ncol=1)
	for (k_ind in (1:dim(netW1)[2])) {
		sum_friendWs <- sum_friendWs + netW1[, k_ind]
	}
  	indivW <- W1
  	A <- rep(0,n)
  	A[which((indivW==1)|(sum_friendWs>cutt_offW))] <- 1 
 	return(A)  	
}	
# Only personal covariates, no interference - DIRECT EFFECT, g(A|W_i)
f.A_nospill <- function(data, ...) plogis(Intercept_A + coeff_A_W*data$W1)
# Only interference effects A_i - INDIRECT EFFECT, g(A|W\W_i)
f.A_spillonly <- function(k, data, ...) {
  	n <- nrow(data)
  	W1 <- data$W1
  	netW1 <- subset(data, select = netW1_1:eval(parse(text=paste("netW1_", k, sep = ""))))
  	sum_friendWs <- matrix(0, nrow=n, ncol=1)
  	# sum_friendW2 <- matrix(0, nrow=n, ncol=1)
	for (k_ind in (1:dim(netW1)[2])) {
		sum_friendWs <- sum_friendWs + netW1[, k_ind] * coeff_A_Wfriend				
	}  	
  	pA <- plogis(Intercept_A + sum_friendWs)
 	return(pA)  	  
}
#---------------------------------------------------------------------------------
# USER DEFINED h_0 and h_star functions
#---------------------------------------------------------------------------------
# f_h_user <- function(k, data, node_l, NetInd_k) {
# 	.f.allCovars <- function(k, NetInd_k, Var, VarNm) {
# 		n <- length(Var) 
# 		NetInd_k <- matrix(NetInd_k, nrow=n, ncol=k)
# 		netVar_names <- NULL
# 		netVar_full <- NULL
# 		d <- matrix(0, nrow=n, ncol = k+1)
# 		d[ , 1] <- Var
# 		d[ , c(2:(k+1))] <- apply(NetInd_k, 2, function(k_indx) {
# 												netVar <- Var[k_indx]
# 												netVar[is.na(netVar)] <- 0
# 												return(netVar)
# 												})
# 		if (k>0) netVar_names <- jpaste(jpaste("net", VarNm, "_"), c(1:k))
# 		Var_names <- c(VarNm, netVar_names)
# 		colnames(d) <- Var_names
# 		return(d)
# 	}
# 	.f.cumprod.matrix <- function(data.indA, data.probA) {
# 		y <- matrix(1, nrow=dim(data.probA)[1], ncol=dim(data.probA)[2])
# 		y[, 1] <- data.probA[,1]^as.integer(data.indA[,1]) *
# 								(1-data.probA[,1])^(1-as.integer(data.indA[,1]))
# 		if (dim(data.probA)[2] > 1) {
# 			for (i in 2:dim(data.probA)[2]) {
# 				y[,i] <- y[,i-1] * (data.probA[,i]^as.integer(data.indA[,i]) * 
# 									(1-data.probA[,i])^(1-as.integer(data.indA[,i])))
# 			}
# 		}
# 		return(round(y[,dim(data.probA)[2]],6))
# 	}	
# 	netAs <- .f.allCovars(k, NetInd_k, data[,node_l$Anode], node_l$Anode)
# 	prob_g0 <- cbind(rep(0.5, nrow(data)), 
# 						apply(NetInd_k, 2, function(k_indx) {
# 											netVar <- rep(0.5, length(k_indx))
# 											netVar[is.na(k_indx)] <- 0
# 											return(netVar)
# 											}))
# 	h_vec_g0 <- .f.cumprod.matrix(netAs, prob_g0)
# 	prob_g_star <- cbind(rep(1, nrow(data)), 
# 						apply(NetInd_k, 2, function(k_indx) {
# 											netVar <- rep(1, length(k_indx))
# 											netVar[is.na(k_indx)] <- 0
# 											return(netVar)
# 											}))

# 	h_vec_g_star <- .f.cumprod.matrix(netAs, prob_g_star)	
# 	return(list(P.hbar.c=list(h_vec=h_vec_g0), P.hbar.star.c=list(h_vec=h_vec_g_star)))
# }
	
#----------------------------------------------------------------------------------
# DATA GENERATING PARAMETERS
# Generate the network using structural equations model
#---------------------------------------------------------------------------------
EC <- 0.35					# Distribution of W, P(W=1)?
NONFIXED_N_FRIENDS <- TRUE	# Simulate # of network friends from uniform (vs fixed)
SYMM_CONN <- FALSE			# Set SYMM_CONN=T to generate symm conn mtx
# coefficients for true Q(Y|A,W)
Intercept_Y <- -2.5 
coeff_Y_W <- 1.5  
# coeff_Y_W <- 2
coeff_Y_Wfriend <- 1.5
# coeff_Y_Wfriend <- 0.7
# old coef. on A in regression for Y
# coeff_Y_A <- 3
coeff_Y_A <- 0.5
coeff_Y_Afriend <- 1.5
 
#Sample 1 community (C_j) with EC and f.g_A_W=A^C
gendata_Cj <- function(C_j = 1, n, k, EC, f.g_name=NULL, f.g_args=NULL) { 
  	#n - number of individuals
  	#k - max size of each individual's network
  	#C_j - community # being sampled (j)
  	#EC - community-level covariate, only influences W_i 
  	#(right now its a probability for binom distr of W_1)
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
  	# SAMPLES SYMMETRIC CONNECTIVITY MATRIX (to the extent possible)  	

  	.f.genConnectMatx_sym <- function(Net_num) {	
	  	I <- big.matrix(n,n, type="short", init=0, shared=FALSE)
		nFriendTot <- array(rep(0,n))
	  	for (index in (1:n)) {
		  		I[index,index] <- 1
		  		FriendSampSet <- setdiff(which(nFriendTot<Net_num), c(0:index)) 
				nFriendSamp <- max(Net_num[index] -	nFriendTot[index], 0)	
				nFriendSamp <- min(nFriendSamp, length(FriendSampSet)) 
				if (nFriendSamp>0) {				
					#To handle case where |FriendSampSet|=1
					if (length(FriendSampSet)==1)  { 
						friends_i <- FriendSampSet	
					} 
					#sample from the possible friend set
					else {	
						friends_i <- sort(sample(FriendSampSet, nFriendSamp))
					}		
					I[index, friends_i] <- 1 
					I[friends_i, index] <- 1
					nFriendTot[index] <- nFriendTot[index] + nFriendSamp
					nFriendTot[friends_i] <- nFriendTot[friends_i] + 1 
				}
		}
		return(I)	
	}
  	# SAMPLES ASYMMETRIC CONNECTIVITY MATRIX (to the extent possible)  	
	.f.genConnectMatx_asym <- function(Net_num) {	
	  	I <- big.matrix(n,n, type="short", init=0, shared=FALSE)
		nFriendTot <- array(rep(0,n))
	  	for (index in (1:n)) {
		  		I[index,index] <- 1
		  		FriendSampSet <- setdiff( c(1:n), index) 
				nFriendSamp <- max(Net_num[index] -	nFriendTot[index], 0)	
				# nFriendSamp <- min(nFriendSamp, length(FriendSampSet)) 
				if (nFriendSamp>0) {
					#To handle case where |FriendSampSet|=1
					if (length(FriendSampSet)==1)  { 
						friends_i <- FriendSampSet	
					} 
					#sample from the possible friend set
					else {	
						friends_i <- sort(sample(FriendSampSet, nFriendSamp))
					}		
					# I[index, friends_i] <- 1 
					I[friends_i, index] <- 1
					nFriendTot[index] <- nFriendTot[index] + nFriendSamp
					# nFriendTot[friends_i] <- nFriendTot[friends_i] + 1 
				}
		}
		return(I)	
	}
  	#Update the # of friends for each individual (given sampled connnectivity matrix)
  	.f.Net_num_update <- function(ConnectMatx) return(colsum(ConnectMatx, c(1:n)) - 1)	
  	# Convert connectivity matx to a vector of friend's ids for each individual i, 
  	.f.Net_vect <- function(ConnectMatx) {	  	
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(which(ConnectMatx[, index]!=0), index)
				netwklist_i	}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}
  	.f.Net_vect_big <- function(ConnectMatx) {
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
				netwklist_i	}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}		
   	.f.Net_vect_bigID <- function(ConnectMatx) {
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
				if (length(netwklist_i)>0) netwklist_i <- paste('I', netwklist_i, sep="")
				return(netwklist_i)}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}		  	
  	.f.mkstrNet <- function(Net) sapply(Net, function(Net_i) paste(Net_i, collapse=" "))	
 	.f_g_wrapper <- function(k, W_netW, fcn_name, ...) {   #wrapper fcn for g(A|W)	
	  	# assign(".f_gAW", get(fcn_name))
	  	args0 <- list(k=k, data=W_netW)
	  	args <- c(args0, ...)
	  	# print("g fcn"); print(fcn_name); print(names(args));
		# A <- rbinom(n, 1, do.call(.f_gAW, args))
		A <- rbinom(n, 1, do.call(fcn_name,args))
		return(A)
	}

   	#get all friends Ws as a list of vectors
  	.f.cA_Pa <-function(W1, Net) sapply(Net, function(netwk) W1[netwk], simplify=F)   
  	.f.redefineCov <- function(Var, VarNm) {
		  	#get all friends Ws as a matrix of dim(n,k) filling unused cols with zeros
		  	.f.netCovar <-function(Covar, Net) sapply(Net, function(netwk) 
		  										c(Covar[netwk],rep(0,k-length(netwk)))) 
			netVar_full <- NULL
			netVar_names <- NULL
			sumnet_name <- NULL
			if (k>0) {
				netVar_full <- .f.netCovar(Var, Net_vect) 
				if (k>1) netVar_full <- t(netVar_full)  
				netVarNm <- paste("net", VarNm, "_", sep="")
				netVar_names <- paste(netVarNm, c(1:k), sep = "") 
			}
			Var_names <- c(VarNm, netVar_names)	
			d <- cbind(Var, netVar_full)			
			colnames(d) <- Var_names
			return(d)   
	}
  	.f.cY_Pa <-function(W1, A, Net) {
  			cY_Pa <- list(Pa_Ws=sapply(Net, function(netwk) W1[netwk], simplify=F),
  			 			Pa_As=sapply(Net, function(netwk) A[netwk], simplify=F)) 
  	}
  	#Define Y, using the network N, W , A, W1_i & A_i
	.f.Y <- function(W1, A, cY_Pa) {  	
		  	f.netY_AorWs <- function(netwk, coeff)  sum(coeff * netwk) 		
		  	sum_friendY_Ws <- sapply(cY_Pa$Pa_Ws, f.netY_AorWs, coeff_Y_Wfriend)
		  	sum_friendY_As <- sapply(cY_Pa$Pa_As, f.netY_AorWs, coeff_Y_Afriend)
		  	indivY_W <- coeff_Y_W * W1
		  	indivY_A <- coeff_Y_A * A
		 	Y <- rbinom(n, 1, plogis(Intercept_Y + sum_friendY_Ws + 
		 								sum_friendY_As + indivY_W + indivY_A))
		 	return(Y)
	}
  #-----------------------------------------------
  #Generating covars
  #-----------------------------------------------
  	print("---------------------------")
  	# print("f.g"); print(f.g_name); print(paste("f.g_args:", f.g_args));
  	W1 <- .f.W1(EC, n)
   	# print(W1)
  	# Set samp=F to assign same # of friends to each person
  	nFriends <- .f.Net_num(samp=NONFIXED_N_FRIENDS, n, k)
  	# print(nFriends)
  	
  	#---------------------------------------------------------
	# Generate symmetric connectivity matrix of friends
	if (SYMM_CONN) ConnectMatx <- .f.genConnectMatx_sym(nFriends) 
	# Generate asymmetric connectivity matrix of friends
		else ConnectMatx <- .f.genConnectMatx_asym(nFriends)  	
  	# print(ConnectMatx[,])  #print the entire connectivity matrix
 	#--------------------------------------------------------- 	
  	nFriends <- .f.Net_num_update(ConnectMatx)  #Update # of actual friends sampled, |N_i|
  	Net_vect <- .f.Net_vect_big(ConnectMatx)    #Get list of vectors of friend's ids for each i
  	Net_vectID <- .f.Net_vect_bigID(ConnectMatx)

  	W_netW <- data.frame(.f.redefineCov(W1, "W1"), nFriends=nFriends)
	# print(head("W_netW")); print(head(W_netW))
  	cA_Pa <- .f.cA_Pa(W1, Net_vect)   #Get c^A - fcn for parents of A_i: (W)
  	# print(cA_Pa)
  	if (is.null(f.g_name)) {
	  	A <- .f_g_wrapper(k, W_netW, f.A) }
	else {
	  	A <- .f_g_wrapper(k, W_netW, f.g_name, f.g_args)
	}

	#convert f.g_args to text format (for output)
  	if (is.null(f.g_args)) {
	  	f.g_args_txt <- "NA" }
	else {
	  		f.g_args_txt <- paste(f.g_args, collapse=" ")
  	}
  	#Get c^Y fcn for parents of Y_i: (A,W)   
  	cY_Pa <- .f.cY_Pa(W1, A, Net_vect)   
  	A_netA <- data.frame(.f.redefineCov(A, "A"))

  	Y <- .f.Y(W1, A, cY_Pa)   
 	#Convert N_i, cA_Pa_i, cY_Pa_i to strings
  	# Net_str <- .f.mkstrNet(Net_vect)   
  	Net_str <- .f.mkstrNet(Net_vectID)

  	netW_str <- .f.mkstrNet(cA_Pa)
  	netA_str <- .f.mkstrNet(cY_Pa$Pa_As)

  	# IDs <- seq(n)  	
	IDs <- paste('I', seq(n), sep='')

  	d <- data.frame(IDs=IDs, C_j, f.g_args_txt, EC, Y, nFriends, W_netW, 
  					A_netA, Net_str, netW_str, netA_str, stringsAsFactors = FALSE)
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
	# print("f.g_list"); print(f.g_list)
	if (is.null(f.g_args_list)) f.g_args_rand<-list(NULL)  
	# d <- do.call("rbind", mapply(gendata_Cj, 1:nC, n_arr, k_arr, EC_rand, 
	  								# f.g_rand, f.g_args_rand, SIMPLIFY=FALSE))
	d <- gendata_Cj(C_j=1, n=n_arr, k=k_arr, EC=EC_rand, f.g_name=f.g_rand, f.g_args=f.g_args_rand)
	d <- d[, (names(d) %in% c('IDs','W1','A','nFriends','Y','Net_str'))]
	return(d)  
  	}

#write.table(d, file = "./netwk_dt.csv", sep = ",")

	# df <-gendata_pop(nC=1, n_arr =10000, k_arr =2, EC_arr=EC, f.g_list=get(f.g0_name), f.g0_args)
	# print(head(df,50))
	# mean(df$Y) # [1] 0.4067
	# df_psi_g.str <-gendata_pop(nC=1, n=10000, k=2, EC, g1_fcns_list_sim, g1_args_list)
	# mean(df_psi_g.str$Y) # [1] 0.5619

#Direct TMLE call
	# df <-gendata_pop(nC=1, n_arr=500, k_arr=2, EC_arr=0.35, f.g_list=get("f.A_x"), f.g_args_list=list(x=0.5))
	# print(head(df,50))
	# mean(df$A)
	# mean(df$Y)
	# df <- df[, (names(df) %in% c('IDs','W1','A','nFriends','Y','Net_str'))]
	# print(head(df, 10))	

			# set.seed(250)
			# n <- narr[1]
			# k <- 2
			# df <-gendata_pop(nC=1, n=200, k=2, EC, f.g0_name, f.g0_args)
			# df <-df[, (names(df) %in% c('IDs','W1','A','nFriends','Y','Net_str'))]
			# # sample_network_k2 <- df[, (names(df) %in% c('IDs','W1','A','nFriends','Y','Net_str'))]
			# # save(sample_network_k2, file="./sample_network_k2.RData")

			# # # get g0 function evaluated
			# f.g0 <- NULL
			# f.g2_args <- NULL
			# if (run_h_true_g0) {
			# 	f.g0 <- get(f.g0_name)
			# 	f.g0_args <- f.g0_args
			# }

			# f.g1 <- get(g1_fcns_list_sim[[1]])
			# args_x <- g1_args_list[[1]]
			# f.g2 <- NULL
			# Qform <- Qform_arr[1]
		  	# gform <- gform_arr[1]
			# tmlenet_out <- tmlenet(data=df, k_max=k, IDnode=IDnode, NETIDnode=NETIDnode,
			# 						Anode=Anode, Wnodes=Wnodes, Ynode=Ynode, nFnode=nFnode,
			# 						Qform=Qform, QDETnode=QDETnode,
			# 						gform=gform, gDETnode=gDETnode, gbound=0.025,
			# 						h_f.g0=f.g0, h_f.g0_args=f.g0_args,
			# 						f.g1.star=f.g1, f.g1_args=args_x,
			# 						f.g2.star=f.g2, f.g2_args=f.g2_args,
			# 						alpha=alpha, h_user_fcn=f_h_user, 
			# 						h_logit_sep_k=h_logit_sep_k)
