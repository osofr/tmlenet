#-----------------------------------------------------------------------------
# ALGORITHM for arbitrary network summary measures of covariates and treatments
#-----------------------------------------------------------------------------
# NOTES.
# *) (W^s,A^s) will be defined by a list of R expressions, captured with substitute() like node() in simcausal

# (I) CREATE 2 data.frames OF SUMMARY MEASURES: sW, sA
      # Approach:
      # Build a data.frame (cYmtx) of all network covariates (A,netA,W,netW), applying .f.allCovars for each W and A.
      # For each expression in (W_i^s[j1], A_i^s[j2]), eval it in the environment of cYmtx, for j1=1,...,s1, j2=1,...,s2.
      # EACH component of W_i^s and A_i^s is a number (column) that is a tranformation of the network data.frame cYmtx
      # e.g., let W=(W1,W2), let W_i^s=(W1, W2, W1*W2 + W1[1]*W2[1] + W1[2]*W1[2] + ..., nFriends) - 4 dimensional summary measure.
      # W1=W1[0] - i's W1, W1[1] - W1 of i's 1st friend, W1[2] - W1 for i's 2nd friends, etc....
      # Run .f.allCovars2 separately for each of the 4 components of W_i^s above, result is a summary measure data.frame of 4 dimensions (columns)
      # Could save memory: run .f.allCovars only for vars that are netVar part of summary measures W_i^s[j], j=1,..s.
      # Could save memory, but longer runtime: recreate cYmtx for each W_i^s[j] and run .f.allCovars only for netVars in W_i^s[j].

# (II) DISCRETIZE ALL CONTINUOUS SUMMARY MEASUREs in sA (do not touch Ws) - SEE BELOW
# A: create NBins dummy indicators
      # Need a good default for discretizing summary measures. Should the number of categories be an input from a user?
      # Bins: 1) Define B_j as equal length intervals on 0-1? or 2) Use quantiles, for an even mass distribution between bins?
      # Continuous - scaled to 0-1
      # Could define indicators as cdf or prob mass function, no difference (see approach 1 & 2 below)!
 
# B (POOLED): create NBins dummy indicators, create a covariate for BIN_j (order of the bin) and pool different BIN_j as P(BIN_j | sW)
      # ...NEED ALGORITHM...

# (IV) FIT P(sA|sW) for summary measures sA[j], j=1,...,s (or P(summaryA[j],j=1,...,s | summaryW))

# (V) FOR GIVEN (sA,sW)=(sa,sw) PREDICT P(sA=sa|sW=sw)?
# Resulting object needs to build a likelihood for (sa,sw) based on discretized models for (sA,sW)...

# (VI) Delayed evaluation for defining subsets 
# (returns a logical vector that will be used for subseting rows of the data)


# ------------------------------------------------------------------------------
# TEST DATASET
# Need a generate a continous summary measure (s.a., rnorm), conditional on covariates
# Test that its correctly estimatated by comparing MSE?
# ------------------------------------------------------------------------------
datatest <- data.frame(
                      W = rbinom(100, 1, prob=.5), 
                      A = rbinom(100, 1, prob=.2),
                      sA = rnorm(100),
                      Y = rbinom(100, 1, prob=.7), 
                      nFriends = rbinom(100, 5, prob=0.2)
                      )
nodes <- list(Anode = "A", Wnode = "W", Ynode = "Y", nFnode = "nFriends")
datatest$sA


# ------------------------------------------------------------------------------
# TEST COMBINATION OF named fun args passed with ... and calling fun with do.call
# ------------------------------------------------------------------------------
insidefun <- function(sA_class, sA_nms, sW_nms, subset, extraparam, ...) {
  print("inside fun: "); print(eval(extraparam))
}
outsidefun <- function(delay.eval = FALSE, ...) {
  # Combine list of params binparams <- list(sA_class, sA_nms, sW_nms, subset) from datnet.sWsA with list(...)
  bin.m.params <- list(sA_class = "test1", sA_nms = "test2", sW_nms = "test3", subset = "test4")

  if (!delay.eval) {
    addl_params <- list(...)  # to evaluate these arguments in the calling environment  
  } else {
    addl_params <- eval(substitute(alist(...))) # to capture the expressions in args and delay their evaluation  
  }

  if (length(addl_params)>0) {
    print("addl_params"); print(addl_params)
  }

  bin.m.params <- append(bin.m.params, addl_params)
  parnames <- names(bin.m.params)
  if (length(bin.m.params) != 0 && (is.null(parnames) || any(parnames==""))) {
    stop("please specify name for each attribute")
  }      
  do.call(insidefun, bin.m.params)
}
c <- 5
d <- 10
outsidefun(delay.eval = FALSE, extraparam = c + d)
outsidefun(delay.eval = TRUE, extraparam = c + d)


# ------------------------------------------------------------------------------
# (II) ALGORITHM FOR DISCRETIZING CONTINUOUS SUMMARY MEASUREs in sA
# ------------------------------------------------------------------------------
  # 1) normalize x to (0,1)
  x <- rnorm(100000)
  # x <- rnorm(100)
  x <- (x - min(x)) / (max(x) - min(x))
  # 2A) turn x into an ordinal (1, 2, 3, ...) based on interval defn 1: bin x by equal length intervals of 0-1:
  intlen <- .1; # generates 10 intervals (intlen <- .2; intlen <- .5;)
  intvec <- seq(from = 0, to = 1, by=intlen)
  system.time(
    i <- data.frame(x = x, xcat = findInterval(x = x, vec = intvec, rightmost.closed = TRUE))
  )
  # for 100K, intlen=.1:
  # user  system elapsed 
  # 0.016   0.001   0.026
  # length(unique(i$xcat)); hist(i$xcat)
  # 2B) turn x into an ordinal (1, 2, 3, ...) based on interval defn 2: bin x by mass (quantiles on 0-1 as probs)
  system.time( {
    quantvec <- quantile(x=x, probs = intvec)
    i["xcatq"] <- findInterval(x = x, vec = quantvec, rightmost.closed = TRUE)    
  })      
  # length(unique(i$xcatq)); hist(i$xcatq)
  # 3) CREATE INDICATOR DUMMIES FOR ALL NON-BINARY SUMMARY MEASUREs in sA (do not touch sW)
  # Approach 1: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
  cats <- sort(unique(i$xcat))
  bvarnmsnew <- "B"%+%"_"%+%cats[-length(cats)]
  system.time({
    dummies2 <- matrix(nrow=nrow(i), ncol=length(cats)-1)
    for(level in cats[-length(cats)]) {
        # dummies2[,level] <- as.integer(i$xcat <= level)
      # rewrite of above (***NOT TESTED***):
      subset_Bj <- i$xcat <= level
      dummies2[subset_Bj, level] <- as.integer(subset_Bj)
      dummies2[!subset_Bj, level] <- gvars$misval
    }
    colnames(dummies2) <- bvarnmsnew
  })
  head(dummies2)
  class(dummies2)
  class(dummies2[,1])
  # for 100K numeric RVs with intlen <- .1;
  #  user  system elapsed 
  # 0.037   0.006   0.054
  # 0.032   0.005   0.038
  # Approach 2 (model.matrix): a bit slower than manual
  # Creates B_j that jumps to 1 and then back to 0 (degenerate) and uses category 1 as ref (droped)
  # First need to make a factor from xcat; adds a column of 1's (removed manually)
  cats <- sort(unique(i$xcat))
  bvarnmsnew <- "B"%+%"_"%+%cats[-1]
  system.time({
    i$xcatf <- factor(i$xcat);
    dummies1 <- model.matrix(~i$xcatf);
    saveattrs <- attributes(dummies1)
    bvarnms <- colnames(dummies1)[-1]
    colnames(dummies1)[-1] <- bvarnmsnew
    # i <- data.frame(i, dummies1[,-1])
  })
  # Approach 2 manual. Creates B_j that jumps to 1 and then back to 0 (degenerate) and exclude the reference category (last)
  cats <- sort(unique(i$xcat))
  bvarnmsnew <- "B"%+%"_"%+%cats[-length(cats)]
  system.time({
    dummies2 <- matrix(nrow=nrow(i), ncol=length(cats)-1)
    for(level in cats[-length(cats)]) {
      dummies2[,level] <- as.integer(i$xcat == level)
    }
    colnames(dummies2) <- bvarnmsnew
  })
