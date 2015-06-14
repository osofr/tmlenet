#------------------------------------
# sim_genk6dat.R
# author: Oleg Sofrygin 
# DATA SIMULATOR FOR NETWORKS WITH K=6 (01/02/2014)
# Source the file to generate a network using these functions
#------------------------------------

#------------------------------------
# SEARCH TRIAL SIMULATION
#------------------------------------			
# Chance of getting infected P(Y=1) is logit linear in 
# no. of infected (W1=1) AND untreated (A=0) friends in the network

# For those who have no infected friends or
# no untreated infected friends P(Y=1)=0

# Thus W's have no direct effect on Y, 
# except for deterministically setting Y=1 if W1=1

# The source of confounding is W2=W_risk, 
# which is predictive of low chance of getting treated, P(A=1), as well as 
# acting on a lot of people (a lot of people have that person in their network)

# THUS, network N is constructed conditional on W2=W_risk

#------------------------------------	
# CONFOUNDING
#------------------------------------	
#1) W1[i] affects the total number of friends i has
#2) W1[i] affects the probability of having [i] selected as someone's friend

# Additionally add individual rate of infection (interceptY)
# Alternative can assume that no one in the network of i is infected & untreated then Risk_i=0
# p_YRisk[sum_friendY_Ws==0] <- 0

#----------------------------------------------------------------------------------
#DEFINE REGRESSION FORMS FORM Q AND g
#----------------------------------------------------------------------------------	
# W1 - inf.RISK
# W2 - (0,1) infected at time t=0     			     
#------------------------------------------------------------------------------ 	
# DEFINE Qform & gform
# Qform Q0 = b0 + b1*I(N(inf & untrt)=0) + b2*N(inf & untrt) + b3*nFriends
# Qform Q0 = b0 + b1*I(sum(netW2 & (1-netA))=0) + b2*N(netW2 & (1-netA)) + 
# b3*nFriends

.f.Qform <- function(k, qform_Vars=NULL, k_miss=1, miss=F) {
	qform_NetVars <- NULL
	if (miss) k <- k_miss
	if (k > 0) { 
		qform_NetVarsW2A <- paste(rep("netW2_",k), c(1: k), "*",
								  rep("(1-netA_",k), c(1: k), ")", sep = "")
		qform_NetVarsW2A <- paste("I(", paste(qform_NetVarsW2A, collapse="+"), ")", sep="")
		qform_NetVars <- paste(qform_NetVarsW2A, collapse="+")
	}
	Qform <- paste("Y ~ ", qform_NetVars, collapse = "")
	return(Qform)
}
# gform  g_0 = b0 + b1*W1 + b2*(sum(netW1))
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

#---------------------------------------------------------------------------------
# DEFINITION OF TREATMENT MECHANISM g(A|W) AND g*(A|W)
#---------------------------------------------------------------------------------
# Actual pop g(A|W) from NPSEM, logistic fcn with a coeff for network Ws and a coeff for individ W_i
# DEFINITION OF g(A_i|W)=P[A=1|cA]:  
# W1 - inf.RISK
# W2 - (0,1) infected at time t=0
# g_0 = b0 + b1*W1 + b2*(sum(netW1))
f.A <- function(k, data, ...) {
  	# coefficients for true g(A|W)
  	# Define coefficient for i individual's W
  	# Intercept_A <- -0.2
  	Intercept_A <- 0.5
  	coeff_A_W <- -0.5
  	coeff_A_Wfriend <- -0.1

  	n <- nrow(data)
  	nFriends <- data$nFriends
  	W1 <- data$W1
  	W2 <- data$W2
	netW1 <- subset(data, select = netW1_1:eval(parse(text=paste("netW1_", k, sep = ""))))
	netW2 <- subset(data, select = netW2_1:eval(parse(text=paste("netW2_", k, sep = ""))))
  	.f_evalA_mtx <- function(W1, W2, netW1, netW2, n, nFriends, ...) {
	  	sum_friendWs <- matrix(0, nrow=n, ncol=1)
		for (k_ind in (1:dim(netW1)[2])) {
			sum_friendWs <- sum_friendWs + netW1[, k_ind] * coeff_A_Wfriend
		}
		return(sum_friendWs)
  	}
  	# P(A) component from individual infection risk (W1)
	indivW <- coeff_A_W * W1
 	sum_friendWs <- .f_evalA_mtx(W1, W2, netW1, netW2, n, nFriends, ...)
  	probA <- plogis(Intercept_A + indivW + sum_friendWs)

  	# set A=0 if W2=0 (not infected at baseline)
  	probA[W2==0] <- 0
 	return(probA)
}
# Set x% of community to A=1 (returns probability P(A=1))
f.A_x <- function(data, x, ...) rep(x, nrow(data))
# Set x% of community to A=1 only among W2=1
  # W1 - inf.RISK
  # W2 - (0,1) infected at time t=0
f.A_xlevelW2_1 <- function(data, x, ...) {
	# print("x"); print(x)
  	n <- nrow(data)
  	W2 <- data$W2
  	pA <- rep(0,n)
	pA[which(W2==1)] <- x
  	return(pA)
}
# Set x% of community to A=1 based on connectivity |N_i|= {low, high}, based on median
f.A_xlevelNi <- function(data, x_Ni_low, x_Ni_high, ...) {	
  	n <- nrow(data)
  	nFriends <- data$nFriends	
    Ni_med = quantile(nFriends, 0.5)
  	x <- rep(0,n)
	x[which(nFriends<= Ni_med)] <- x_Ni_low
	x[which(nFriends> Ni_med)] <- x_Ni_high
  	return(x)  
} 	 
# Deterministically set A=1 based on cutt-off value for Ws or if W_i = 1
f.A_cutt_offW <- function(k, data, cutt_offW, ...) {	
  	n <- nrow(data)
  	W1 <- data$W1
  	netW1 <- subset(data, select = netW1_1:eval(parse(text=paste("netW1_", k, 
  																	sep = ""))))
  	sum_friendWs <- matrix(0, nrow=n, ncol=1)
	for (k_ind in (1:dim(netW1)[2])) {
		sum_friendWs <- sum_friendWs + netW1[, k_ind]
	}
  	indivW <- W1
  	A <- rep(0,n)
  	A[which((indivW==1)|(sum_friendWs>cutt_offW))] <- 1 
 	return(A)  	
}
 
#---------------------------------------------------------------------------------
# Generate the network population using structural equations model
#---------------------------------------------------------------------------------
EC <- 0.35
NONFIXED_N_FRIENDS <- TRUE 	# Simulate # of network friends from uniform (vs fixed)
SYMM_CONN <- FALSE 	# Set SYMM_CONN=T to generate symmetric connectivity mtx
# GETS INCIDENCE [1] 0.2278, k=10 & [1] 0.192, k=6 
  # Intercept_Y <- -3.7
  # # Risk Y=1 decrement if no one in i's network is infectous
  # coeff_Y_AW_0 <- -1
  # # Risk Y=1 increment for every additional infectous partner
  # # coeff_Y_AW <- 0.9
  # coeff_Y_AW <- 2
# GETS INCIDENCE: 0.0394, k=6 
Intercept_Y <- -4.7
  # Risk Y=1 decrement if no one in i's network is infectous
coeff_Y_AW_0 <- -1
  # Risk Y=1 increment for every additional infectous partner
  # coeff_Y_AW <- 0.9
# coeff_Y_AW <- 1.1
coeff_Y_AW <- 2

#Sample 1 community (C_j) with EC and f.g_A_W=A^C
gendata_Cj <- function(C_j = 1, n, k, EC, f.g_name=NULL, f.g_args=NULL) { 
  	#n - number of individuals
  	#k - max size of each individual's network
  	#C_j - community # being sampled (j)
  	#EC - community-level covariate, only influences W_i
  	#f.g_name - fcn for community intervention on A^C, i.e. g(A|W)
  	#f.g_args - additional args to be passed to f.g_name (as a list)

 	#----------------------------------------------------------------
  	# Defining structural equations 
  	#----------------------------------------------------------------
  	
  	# W1 - categorical or continuous confounder (5 categories, 0-4)
  	# W1 - risk of future infection
  	nW1cat <- 6
  	.f.W1 <- function(n) rbinom(n, 5, prob=c(0.4, 0.5, 0.7, 0.4))
	
  	# W2 - binary infection status at t=0, positively correlated with W1
  	.f.W2 <- function(Cj_prob, W1, n) {
  		# prob_W2 <- seq(0.25, 0.5, by=0.3/nW1cat)
  		prob_W2 <- seq(0.45, 0.8, by=0.3/nW1cat)
  		return(sapply(c(1:n), function(i) rbinom(1, 1, prob=prob_W2[W1[i]+1])))
  	}
	# Total number of friends for each i, influenced by W1 - inf. risk
  	.f.Net_num <- function(W1, samp=FALSE, n, k) {
		# network confounding through W2
		k_arr <-c(1:k)
		pN_0 <- 0.02
		#1) W1[i] affects the total number of friends i has (prob goes up)
		# W1=0 probabilities of |F_i|
		prob_Ni_W1_0 <- c(pN_0, plogis(-3 - 0 - k_arr/2))		
		# W1=1 probabilities of |F_i|
		prob_Ni_W1_1 <- c(pN_0, plogis(-1.5 - 0 - k_arr/3))		
		# W1=2 probabilities of |F_i|
		prob_Ni_W1_2 <- c(pN_0, pnorm(-2*abs(2 - k_arr) / 5))	
		# W1=3 probabilities of |F_i|
		prob_Ni_W1_3 <- c(pN_0, pnorm(-2*abs(3 - k_arr) / 5)	)	
		# W1=4 probabilities of |F_i|
		prob_Ni_W1_4 <- c(pN_0, plogis(-4 + 2*(k_arr-2)))		
		# W1=5 probabilities of |F_i|
		prob_Ni_W1_5 <- c(pN_0, plogis(-4 + 2*(k_arr-3)))
		
		prob_Ni <- list(prob_Ni_W1_0, prob_Ni_W1_1,
						prob_Ni_W1_2, prob_Ni_W1_3, 
						prob_Ni_W1_4, prob_Ni_W1_5)
		prob_Ni <- lapply(prob_Ni, function(x) x/sum(x))

		Net_num <- sapply(c(1:n), function(i) sample(0:k, 1, replace=T, 
														prob=prob_Ni[[W1[i]+1]]))
  		return(Net_num)
  	}

  	# W1-based probs of i being selected as someone's friend
	W1cat_arr <- c(1:nW1cat)/2
  	prob_f <- plogis(-4.5 + 2.5*W1cat_arr) / sum(plogis(-4.5 + 2.5*W1cat_arr))
  	
  	# Sample connectivity matrix (0,1), 
  	# randomly selecting from the list of available individuals 
  	# so that we never exceed |N_i| for each i 
  	# may result in actual # of friends being < |N_i| for some
	.f.genConnectMatx_asym_biased <- function(W1, Net_num) {	
	  	I <- big.matrix(n,n, type="short", init=0, shared=FALSE)
		nFriendTot <- array(rep(0,n))
	  	for (index in (1:n)) {
	  		I[index,index] <- 1
			#set of possible friends to sample, anyone but itself		  		
	  		FriendSampSet <- setdiff( c(1:n), index) 
			#check i's network is not already filled to max	  		
			nFriendSamp <- max(Net_num[index] -	nFriendTot[index], 0)	
			if (nFriendSamp>0) {
				#To handle case where |FriendSampSet|=1
				if (length(FriendSampSet)==1)  { 
					friends_i <- FriendSampSet	
				} 
				#sample from the possible friend set, with prob based on W2=W_risk
				else {	
				#2) W1[i] affects the probability of having [i] selected as someone's friend
				# W1-based probs of i being selected
					friends_i <- sort(sample(FriendSampSet, nFriendSamp, prob=prob_f[W1[FriendSampSet]+1]))
				}		
				I[friends_i, index] <- 1
				nFriendTot[index] <- nFriendTot[index] + nFriendSamp
			}
		}
		return(I)			
	}
  	#Update the # of friends for each individual (given sampled connnectivity matrix)
  	.f.Net_num_update <- function(ConnectMatx) return(colsum(ConnectMatx, c(1:n)) - 1)	
  	#Convert connectivity matx to a lists of vector of friend's ids
  	.f.Net_vect <- function(ConnectMatx) {	  	
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(which(ConnectMatx[, index]!=0), index)
				netwklist_i	}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}
	#Same, but using mwhich() instead of which() (faster for n>30,000)	
  	.f.Net_vect_big <- function(ConnectMatx) {
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
				netwklist_i	}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}		
  	.f.Net_vect_bigID <- function(ConnectMatx) {
		f.netwklist_i <- function(index)  {
				netwklist_i <- setdiff(mwhich(ConnectMatx, index, 0, 'neq'), index)
				if (length(netwklist_i)>0) netwklist_i <- paste('I', netwklist_i,
				 															sep="")
				return(netwklist_i)}
		sapply(1:n, f.netwklist_i, simplify=F)	
	}		  	
  	.f.mkstrNet <- function(Net) sapply(Net, function(Net_i) paste(Net_i, 
  																collapse=' '))
  	.f.mkstrNetWi <- function(Net,i) sapply(Net, function(Net_i) paste(Net_i[[i]], 
  																collapse=" "))
  	
  	.f_g_wrapper <- function(k, W_netW, fcn_name, ...) {   #wrapper fcn for g(A|W)	
	  	# assign(".f_gAW", get(fcn_name))
	  	args0 <- list(k=k, data=W_netW)
	  	args <- c(args0, ...)
		A <- rbinom(n, 1, do.call(fcn_name, args))
	}
  	# get all friends Ws as a list of vectors
  	.f.cA_Pa <-function(W1, W2, Net) sapply(Net, function(netwk) 
  											list(W1=W1[netwk], W2=W2[netwk]),
  																	 simplify=F)
  	# get network A's & W's as a matrix
  	.f.redefineCov <- function(Var, VarNm) {
	  	#get all friends Ws as a matrix of dim(n,k) filling unused cols with zeros
	  	.f.netCovar <-function(Covar, Net) sapply(Net, function(netwk) 
	  											c(Covar[netwk],rep(0,
	  														k-length(netwk))))
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
  	# Calculate c^Y(Pa(Y_i)) - a function into R^(k+1), 
  	# doesn't depend on n or i and is permutation invariant
  	# for each i, get the W's and A's in the network => cY is a vector of (W,A),
  	# not including individual (W_i,A_i)
  	.f.cY_Pa <-function(W1, W2, A, Net) {
		cY_Pa <- list(Pa_Ws=sapply(Net, function(netwk) list(W1=W1[netwk], 
																W2=W2[netwk]),
  			 						simplify=F),
  						Pa_As=sapply(Net, function(netwk) A[netwk], simplify=F))
  	}
  	#Define Y, using the network size |N|, netW's , netA, W_i's, A_i
  	# Q0 = b0 + b1*I(N(inf & untrt)=0) + b2*N(inf & untrt) + b3*N_friends
	.f.Y <- function(W1, W2, A, cY_Pa) {  
	  	# GET (total # of friends who are infected AND untreated):
	  	f.netY_AWs <- function(i, coeff)  sum(coeff * cY_Pa$Pa_Ws[[i]]$W2 * (1-cY_Pa$Pa_As[[i]]))
	  	sum_friendY_Ws <- sapply(c(1:n), f.netY_AWs, coeff_Y_AW)

	  	# # N untreated in i's network
	  	# f.netY_As <- function(i)  sum(1-cY_Pa$Pa_As[[i]])
	  	# sum_friendA0 <- sapply(c(1:n), f.netY_As)
	  	# # N infected in i's network		
	  	# f.netY_W2s <- function(i)  sum(cY_Pa$Pa_Ws[[i]]$W2)
	  	# sum_friendW21 <- sapply(c(1:n), f.netY_W2s)
	  	
	  	p_YRisk <- plogis(Intercept_Y + coeff_Y_AW_0*(sum_friendY_Ws==0) + sum_friendY_Ws)
	 	Y <- rbinom(n, 1, p_YRisk)

		# Set outcome Y=0 for those who have W2=1 (infected at baseline)
	 	Y[W2==1] <- 0
	 	return(Y) 
	}
  #-----------------------------------------------
  #Generating covars
  #-----------------------------------------------  
  	print("---------------------------")
  	# print("f.g:"); print(f.g_name); 
  	print(paste("f.g_args:", f.g_args));

  #-----------------------------------------------    	
  	W1 <- .f.W1(n)  	# categorical infection risk
  	W2 <- .f.W2(EC, W1, n) # infected at baseline
	# Generate # of friends, |N_i|
  	nFriends <- .f.Net_num(W1, samp=NONFIXED_N_FRIENDS, n, k)
  	#---------------------------------------------------------
	# Generate symmetric connectivity matrix of friends, by infection risk W1
	if (SYMM_CONN) ConnectMatx <- .f.genConnectMatx_sym(nFriends) else 
					ConnectMatx <-  .f.genConnectMatx_asym_biased(W1, nFriends)  	
	# Generate asymmetric connectivity matrix of friends
  	# print(ConnectMatx[,])  #print the entire connectivity matrix
 	#--------------------------------------------------------- 	
 	#Update # of actual friends sampled, |N_i|
  	nFriends <- .f.Net_num_update(ConnectMatx) 
  	 #Get list of vectors of friend's ids for each i 
  	Net_vect <- .f.Net_vect_big(ConnectMatx)   
  	Net_vectID <- .f.Net_vect_bigID(ConnectMatx)

  	W1_netW <- data.frame(.f.redefineCov(W1, "W1"))
   	W2_netW <- data.frame(.f.redefineCov(W2, "W2"))
   	W_netW <- cbind(W1_netW, W2_netW)
	# print(head("W_netW")); print(head(W_netW))
  	cA_Pa <- .f.cA_Pa(W1, W2, Net_vect)   #Get c^A - fcn for parents of A_i: (W)
  	
  	if (is.null(f.g_name)) {  
	  	A <- .f_g_wrapper(k, W_netW, "f.A") }
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
  	cY_Pa <- .f.cY_Pa(W1, W2, A, Net_vect)   
  	A_netA <- data.frame(.f.redefineCov(A, "A"))

  	Y <- .f.Y(W1, W2, A, cY_Pa)   
 	#Convert N_i, cA_Pa_i, cY_Pa_i to strings
  	# Net_str <- .f.mkstrNet(Net_vect)
  	Net_str <- .f.mkstrNet(Net_vectID)
  	netW1_str <- .f.mkstrNetWi(cA_Pa, 1)
  	netW2_str <- .f.mkstrNetWi(cA_Pa, 2)  	
  	netA_str <- .f.mkstrNet(cY_Pa$Pa_As)
 
 	# Flag for no risk of infection (no infected untreated partners)
  	NoRskF <- sapply(c(1:n), function(i) sum(cY_Pa$Pa_Ws[[i]]$W2 * (1-cY_Pa$Pa_As[[i]]))==0)
	
	IDs <- paste('I', seq(n), sep='')
 	# Add N(infected & untreated) and net1(inf&untrt), ..., netk(inf&untrt) 	
  	d <- data.frame(IDs=IDs, C_j, f.g_args_txt, EC, Y, nFriends, W_netW, 
  					A_netA, NoRskF=NoRskF, Net_str, netW1_str, netW2_str,
  					netA_str, stringsAsFactors = FALSE)
	#***************************   	
	# 1) Set A=0 when W2=0 (exclude W2=0 from modelling g_0)
	# 2) Set Y=1 when W2=1 (exclude W2=1 from modelling Q_0)   	
	#***************************
	d$g_deterministic <- (d$W2==0)
	d$Q_deterministic <- (d$W2==1)
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
	d <- d[, (names(d) %in% c('IDs','W1','W2','A','nFriends','Y','Net_str', 
										'g_deterministic', 'Q_deterministic'))]
	return(d)
}


#------------------------------------------------------------------------------ 
# TO GENERATE A DATASET:
#------------------------------------------------------------------------------ 
# library(bigmemory)
# library(biganalytics)
# options(bigmemory.typecast.warning=FALSE)

# n <- 5000
# k <- 6
# f.g0_name = "f.A"
# f.g0_args = list(NULL)
# df <- gendata_pop(nC=1, n, k, EC, f.g0_name, f.g0_args)
# head(df)

#------------------------------------------------------------------------------ 
# 1) PLAY WITH GENERATING SOME DATA
# 1 YR INCIDENCE OF HIV (2011):
# > 375/100000 
# [1] 0.00375
# 0.00375*100=0.375
# TOTAL: 0.4% percent or 0.004 probability in a year
# 2) CHECK THERE IS CONFOUNDING PRESENT!!!!! -check
#------------------------------------------------------------------------------ 
# options(echo=TRUE)
# .f.allCovars <- function(k, NetInd_k, Var, VarNm) {
# 	n <- length(Var) 
# 	NetInd_k <- matrix(NetInd_k, nrow=n, ncol=k)
# 	netVar_names <- NULL
# 	netVar_full <- NULL
# 	d <- matrix(0, nrow=n, ncol = k+1)
# 	d[ , 1] <- Var
# 	d[ , c(2:(k+1))] <- apply(NetInd_k, 2, function(k_indx) {
#                     											netVar <- Var[k_indx]
#                     											netVar[is.na(netVar)] <- 0
#                     											return(netVar)
#                     											})
# 	if (k>0) netVar_names <- jpaste(jpaste("net", VarNm, "_"), c(1:k))
# 	Var_names <- c(VarNm, netVar_names)
# 	colnames(d) <- Var_names
# 	return(d)
# }
# .f.cumprod.matrix <- function(data.indA, data.probA) {
# 	y <- matrix(1, nrow=dim(data.probA)[1], ncol=dim(data.probA)[2])
# 	y[, 1] <- data.probA[,1]^as.integer(data.indA[,1]) *
# 							(1-data.probA[,1])^(1-as.integer(data.indA[,1]))
# 	if (dim(data.probA)[2] > 1) {
# 		for (i in 2:dim(data.probA)[2]) {
# 			y[,i] <- y[,i-1] * (data.probA[,i]^as.integer(data.indA[,i]) * 
# 								(1-data.probA[,i])^(1-as.integer(data.indA[,i])))
# 		}
# 	}
# 	return(round(y[,dim(data.probA)[2]],6))
# }
# .f.gen.probA_N <- function(df, deterministic, m.gN) {
#     g_N <- predict(m.gN, newdata=df, type="response")
#   	g_N[deterministic] <- 0
#   	return(g_N)
# }
# .f.gen.probA.star <- function(k, df_AllW, fcn_name, f_args=NULL) {	 	
# 	.f_g_wrapper <- function(k, df_AllW, fcn_name, ...) {	
#   		args0 <- list(k=k, data=df_AllW)
#   		args <- c(args0, ...)
#   		# print("f.gen.probA.star fcn_name"); print(fcn_name)
# 		do.call(fcn_name, args)
# 	}
#   	.f_g_wrapper(k, df_AllW, fcn_name, f_args)
# }

# iptw_est <- function(k, data, m.gN, f.g.star, f.g_args, NetInd_k) {	
# 	family="binomial"
# 	n <- nrow(data)
# 	W1names <- c("W1", "netW1_1", "netW1_2", "netW1_3", "netW1_4", "netW1_5", "netW1_6")
# 	W2names <- c("W2", "netW2_1", "netW2_2", "netW2_3", "netW2_4", "netW2_5", "netW2_6")
# 	netW <- data[,c(W1names, W2names)]

# 	Anames <- c("A", "netA_1", "netA_2", "netA_3", "netA_4", "netA_5", "netA_6")
# 	netA <- data[,Anames]

# 	df_AllWs <- netW
# 	cY.mtx <- cbind(netW, netA, subset(data, select="nFriends"))
# 	# Get A's
# 	indA <- netA
#   	determ.g <- data$g_deterministic
# 	# predict g*(A=1|W)
# 	print(f.g_args)
# 	# print(head(cY.mtx))
# 	pA_gstar <- .f.gen.probA.star(k, cY.mtx, f.g.star, f.g_args)
# 	# print(head(pA_gstar, 10))
# 	netpA_gstar <- .f.allCovars(k, NetInd_k, pA_gstar, "A")
# 	# print(head(netpA_gstar, 10))

# 	# calculate likelihoods P_*(A=a|W)
# 	gstar_A <- .f.cumprod.matrix(indA, netpA_gstar)

# 	# predict g0_N(A=1|W)
# 	pA_g0N <- .f.gen.probA_N(cY.mtx, determ.g, m.gN)
# 	netpA_g0 <- .f.allCovars(k, NetInd_k, pA_g0N, "A")
# 	# print("net P_0(A=1)"); print(head(netpA_g0))
# 	# calculate likelihoods P_0(A=a|W)
# 	g0_A <- .f.cumprod.matrix(indA, netpA_g0)
# 	# print("likelihood g(A=a)"); print(head(g0_A))
#   	ipweights <- gstar_A / g0_A
# 	return(ipweights)
# }

# #------------------------------------------------------------------------------
# # Netwk ids strings to lists of vectors
# .f.mkvecNetID <- function(Net_str) lapply(Net_str, function(Net_str_i)
#     										unlist(strsplit(Net_str_i, ' ',
#     										fixed=TRUE)))
# # Netwk ids strings to arr of ID vectors (filled up with trailing NA's)
# .f.mkarrNetID <- function(Net_str) t(sapply(Net_str, function(Net_str_i) {
#     								netwk <- unlist(strsplit(Net_str_i, ' ',
#     								fixed=TRUE))
#     								return(c(netwk, rep(NA,k-length(netwk))))
#     								} ))
# .f_NetInd_k <- function(df) {
# 	# get list of network ID vectors from network strings for each i
# 	NetIDVec <- .f.mkvecNetID(df[,NETIDnode])
# 	# convert into list of network row #s
# 	NetVec <- lapply(NetIDVec, function(NetID) as.numeric(sapply(NetID, function(x) which(x==as.vector(df[,IDnode])))))
# 	# get array of network IDs
# 	NetID_k <- .f.mkarrNetID(df[,NETIDnode])
# 	# make NetID_k into array of network indices (row #s)
# 	NetInd_k <- apply(NetID_k, 2, function(k_ID) {
# 	                sapply(as.vector(k_ID), function(x) {
# 	                  if (is.na(x)) {
# 	                    NA
# 	                  } else { 
# 	                    which(x==as.vector(df[,IDnode]))
# 	                  }
# 	                 })
# 	              })
# 	head(NetInd_k, 50)	
# 	return(NetInd_k)
# }

# 			# n <- n_arr[1]
# 			# source("tmlenet.R")
# 			n <- 10000
# 			k <- 6
# 			f.g0_name = "f.A"
# 			f.g0_args = list(NULL)
# 			df <- gendata_pop(nC=1, n, k, EC, f.g0_name, f.g0_args)
# 			head(df)
# 			# head(cbind(df$A, df$Y), 100)
# 			# probA <- f.A(6, df)
# 			# tapply(df$Y, df$A, mean)
# 			mean(df$nFriends)
# 			tapply(df$nFriends, df$W1, mean)
# 			table(df$W1)/n
# 			mean(df$A)
# 			mean(df$Y)
# 			mean(df$W2)
# 			tapply(df$Y, df$W2, mean)
# 			tapply(df$A, df$W1, mean)
# 			tapply(probA, df$W1, mean)
# 			tapply(df$Y, df$W1, mean)
# 			fit_glm0 <- glm("A ~  W1 + netW1_1 + netW1_2 + netW1_3 + netW1_4 + netW1_5 + netW1_6 + nFriends",
# 							family="binomial", data=df)
# 			fit_glm1 <- glm("A~W1", family="binomial", data=df)
# 			fit_glm3 <- glm("A~W2", family="binomial", data=df)
# 			fit_glm2 <- glm("A~nFriends", family="binomial", data=df)
# 			fit_glm4 <- glm("A~netW2_2", family="binomial", data=df)
			
# 			# compare true prob(A) with predicted based on the glm model
# 			fit_glm1 <- glm("A~W1", family="binomial", data=df)
# 			probA_pred <- predict(fit_glm1, type = "response")
# 			probA_0 <- f.A(k, df)
# 			tapply(probA_pred, df$W1, mean)
# 			tapply(probA_0, df$W1, mean)

# 			tapply(probA_pred, df$netW1_6, mean)
# 			tapply(probA_0, df$netW1_6, mean)
# 			tapply(df$Y, df$netW1_1, mean)



# NETIDnode <- 'Net_str'
# IDnode <- 'IDs'
# determ.Q <- df[,"Q_deterministic"]
# NetInd_k <- .f_NetInd_k(df)

# f.gstar = "f.A_xlevelW2_1"
# f.gstar_args = list(x=0)
# df_true <- gendata_pop(nC=1, n, k, EC, f.gstar, f.gstar_args)

# p_wt_ests <- iptw_est(k=6, data=df, m.gN=fit_glm2, f.g.star=f.gstar, f.g_args=f.gstar_args, NetInd_k=NetInd_k) 
# # head(p_wt_ests)
# Y_IPTW_net <- rep(0, n)
# Y_IPTW_net[!determ.Q] <- df$Y[!determ.Q] * p_wt_ests[!determ.Q]
# # % total infected
# mean(df_true$Y)
# mean(Y_IPTW_net)
# # % of newly infected:
# tapply(df_true$Y, df_true$W2, mean)
# tapply(Y_IPTW_net, df$W2, mean)



			# df <- df[, (names(df) %in% c('IDs','W1','W2','A','nFriends','Y','Net_str',
			# 'g_deterministic', 'Q_deterministic'))]
			# head(df)
			# mean(df$nFriends)

			# # get g0 function evaluated
			# f.g0 <- NULL
			# if (run_h_true_g0) {
			# 	f.g0 <- get(f.g0_name)
			# 	f.g0_args <- f.g0_args
			# }
			# f.g1 <- get(f.g1star_name)
			# f.g2 <- NULL
			# if (ATE_param) {
			# 	f.g2 <- get(f.g2star_name)
			# }
			# Qform <- Qform_arr[1]
  			# gform <- gform_arr[1]

		 	#  	# Get estimators for each x% of intervention g*
		 	#  	args_x <- fargs_vec[1]
		 	#  	psi_0 <- 0
		 	#  	h_logit_sep_k <- TRUE

			# tmlenet_out <- tmlenet(data=df, k_max=k, IDnode=IDnode, NETIDnode=NETIDnode,
			# 						Anode=Anode, Wnodes=Wnodes, Ynode=Ynode, nFnode=nFnode,
			# 						Qform=Qform, QDETnode=QDETnode,
			# 						gform=gform, gDETnode=gDETnode, gbound=0.025,
			# 						h_f.g0=f.g0, h_f.g0_args=f.g0_args,
			# 						f.g1.star=f.g1, f.g1_args=args_x,
			# 						f.g2.star=f.g2, f.g2_args=f.g2_args,
			# 						alpha=alpha, h_user_fcn=f_h_user, 
			# 						h_logit_sep_k=h_logit_sep_k)