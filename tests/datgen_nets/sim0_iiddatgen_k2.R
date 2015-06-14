#------------------------------------
# sim0_iiddatgen_k2
# author: Oleg Sofrygin 
# iid DATA (W,A,Y) with NETWORK K=2 (network has no effect on anything)
# Source the file to generate a network using these functions
#------------------------------------

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
	#---------------------------------------------------------------------------------
  	n <- nrow(data)
  	W1 <- data$W1
  	# P(A) component from individual infection risk (W1)
	indivW <- coeff_A_W * W1
  	probA <- plogis(Intercept_A + indivW)
 	return(probA)   	
}
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
#----------------------------------------------------------------------------------
# DATA GENERATING PARAMETERS
# Generate the network using structural equations model
#---------------------------------------------------------------------------------
NONFIXED_N_FRIENDS <- TRUE	# Simulate # of network friends from uniform (vs fixed)
SYMM_CONN <- FALSE			# Set SYMM_CONN=T to generate symm conn mtx
# coefficients for true Q(Y|A,W)
Intercept_Y <- -2.5 
coeff_Y_W <- 1.5  
coeff_Y_Wfriend <- 0
# coeff_Y_A <- 0.5
coeff_Y_A <- 0.9
coeff_Y_Afriend <- 0

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
  	W1 <- .f.W1(EC, n)
  	nFriends <- .f.Net_num(samp=NONFIXED_N_FRIENDS, n, k)  	# Set samp=F to assign same # of friends to each person
	if (SYMM_CONN) ConnectMatx <- .f.genConnectMatx_sym(nFriends)  # Generate symmetric connectivity matrix of friends
		else ConnectMatx <- .f.genConnectMatx_asym(nFriends) # Generate asymmetric connectivity matrix of friends
  	nFriends <- .f.Net_num_update(ConnectMatx)  #Update # of actual friends sampled, |N_i|
  	Net_vect <- .f.Net_vect_big(ConnectMatx)    #Get list of vectors of friend's ids for each i
  	Net_vectID <- .f.Net_vect_bigID(ConnectMatx)
  	W_netW <- data.frame(.f.redefineCov(W1, "W1"), nFriends=nFriends)
  	cA_Pa <- .f.cA_Pa(W1, Net_vect)   #Get c^A - fcn for parents of A_i: (W)
  	if (is.null(f.g_name)) {
	  	A <- .f_g_wrapper(k, W_netW, f.A) }
	else {
	  	A <- .f_g_wrapper(k, W_netW, f.g_name, f.g_args)
	}
	#convert f.g_args to text format (for output):
  	if (is.null(f.g_args)) {
	  	f.g_args_txt <- "NA" }
	else {
	  		f.g_args_txt <- paste(f.g_args, collapse=" ")
  	}
  	
  	cY_Pa <- .f.cY_Pa(W1, A, Net_vect)   #Get c^Y fcn for parents of Y_i: (A,W)   
  	A_netA <- data.frame(.f.redefineCov(A, "A"))

  	Y <- .f.Y(W1, A, cY_Pa)   
  	Net_str <- .f.mkstrNet(Net_vectID) #Convert N_i, cA_Pa_i, cY_Pa_i to strings
  	netW_str <- .f.mkstrNet(cA_Pa)
  	netA_str <- .f.mkstrNet(cY_Pa$Pa_As)
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
	if (is.null(f.g_args_list)) f.g_args_rand<-list(NULL)  
	d <- gendata_Cj(C_j=1, n=n_arr, k=k_arr, EC=EC_rand, f.g_name=f.g_rand, f.g_args=f.g_args_rand)
	d <- d[, (names(d) %in% c('IDs','W1','A','nFriends','Y','Net_str'))]
	return(d)  
  	}