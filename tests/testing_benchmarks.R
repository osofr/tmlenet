library(pryr)

#---------------------------------------------------------------------
### TESTING BinarySummary and DatSummary clases ###
#---------------------------------------------------------------------
# simulate data, create instances of classes and fit models:
nsamp <- 100
outvarnms <- c("A", "A_net1", "A_net2")
predvarnms <- "W"%+%(1:10)
# make up some data and test the functions:
gen_testdat <- function(predvarnms) {
  testdat <- data.frame(
      A = rbinom(n = nsamp, size = 1, prob = 0.3),
      A_net1 = rbinom(n = nsamp, size = 1, prob = 0.6),
      A_net2 = rbinom(n = nsamp, size = 1, prob = 0.1),
      lapply(predvarnms, function(name) rbinom(n = nsamp, size = 1, prob = 0.3)))
  names(testdat) <- c(outvarnms, predvarnms)
  testdat
}
testdat <- gen_testdat(predvarnms)
ncol(testdat)
dim(testdat)
names(testdat)

summeas <- Summaries$new(Kmax = 2, sA_nms = outvarnms, sW_nms = predvarnms)
summeas$regsnms
summeas$subset_regs
# summeas$deterministic(testdat)
summeas$fit(data = testdat)
summeas$getfits()

# test predict in BinarySummary
lapply(summeas$getfitted.pbins(), function(obj) obj$is.fitted)
lapply(summeas$getfitted.pbins(), function(obj) obj$getprobA1())
lapply(summeas$getfitted.pbins(), function(obj) obj$predictAeqa(testdat$A))
lapply(summeas$getfitted.pbins(), function(obj) obj$getprobAeqa())

# test predictAeqa in BinarySummary
summeas$predict(newdata = testdat)
# test predict and predictAeqa in Summaries
summeas$predictAeqa(indA = as.matrix(testdat[,outvarnms]))
summeas$getcumprodAeqa()

#---------------------------------------------------------------------
### TESTING BinarySummary and DatSummary clases ###
#---------------------------------------------------------------------
# simulate data, create instances of classes and fit models:
nsamp <- 100000
outvarnm <- "A"
predvarnms <- "W"%+%(1:10)
# make up some data and test the functions:
gen_testdat <- function(predvarnms) {
  testdat <- data.frame(
      A = rbinom(n = nsamp, size = 1, prob = 0.3), 
      lapply(predvarnms, function(name) rbinom(n = nsamp, size = 1, prob = 0.3)))
  names(testdat) <- c(outvarnm, predvarnms)
  testdat
}
testdat <- gen_testdat(predvarnms)
ncol(testdat)
dim(testdat)
names(testdat)

# use glm.fit:
datbin_glm <- DatSummary$new(data = testdat, outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp))
class(datbin_glm)
datbin_glm$show()
pbin1 <- BinarySummary$new(outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp))
pbin1$show()
pbin1$is.fitted
pbin1$fit(data_obj = datbin_glm)

# use speedglm.wfit:
datbin_sglm <- DatSummary$new(data = testdat, outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp), glm = FALSE)
class(datbin_sglm)
pbin2 <- BinarySummary$new(outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp))
pbin2$show()
pbin2$is.fitted
pbin2$fit(data_obj = datbin_sglm)
pbin2$getfit()
all.equal(pbin1$getfit()$coef, pbin2$getfit()$coef)


fun.pbinglm <- function() {
  testdat <- gen_testdat(predvarnms)
  datbin_glm <- DatSummary$new(data = testdat, outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp))
  pbin <- BinarySummary$new(outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp))
  pbin$fit(data_obj = datbin_glm)
}
fun.pbinspeedglm <- function() {
  testdat <- gen_testdat(predvarnms)
  datbin_sglm <- DatSummary$new(data = testdat, outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp), glm = FALSE)
  pbin <- BinarySummary$new(outvar_chr = outvarnm, predvars_chr = predvarnms, subsetidx = rep(TRUE,nsamp))
  pbin$fit(data_obj = datbin_sglm)
}

library(microbenchmark)
nsamp <- 100000
speed <- microbenchmark(times=50,
  fun.pbinglm(),
  fun.pbinspeedglm()
)
speed
(speed <- mb_summary(speed))



library(Matrix)
library(speedglm)
regform <- as.formula("A~"%+%paste(predvarnms, collapse="+"))

# vs. 1a/b cbinding from data.frame
fdesmat1 <- function() { as.matrix(cbind(Intercept = 1, testdat[,predvarnms])) }
# vs. 3 pre-allocating design mat w/ loop
fdesmat2 <- function() {
  npred <- length(predvarnms)
  desmat3 <- matrix(1, nrow = nrow(testdat), ncol = npred+1)
  for (coli in seq(predvarnms)) {
    desmat3[,coli+1] <- testdat[,predvarnms[coli]]
  }
  desmat3
}

library(microbenchmark)
speed <- microbenchmark(
  fdesmat1(),
  fdesmat2()
)
mb_summary <- function(x) {
  res <- summary(x, unit="us")
  data.frame(name = res$expr, median = res$median)
}
speed
(speed <- mb_summary(speed))
# Unit: milliseconds
#        expr       min        lq      mean    median        uq       max neval cld
#  fdesmat1()  3.110395  3.602966  6.878445  6.361819  8.004027  17.83923   100  a 
#  fdesmat2() 24.217721 28.027727 31.788855 29.450551 33.214587 131.82055   100   b
# > (speed <- mb_summary(speed))
#         name    median
# 1 fdesmat1()  6361.818
# 2 fdesmat2() 29450.551

system.time(m.glm <- glm(formula=regform, data=testdat, family = binomial()))
m.glm$coef
  # user  system elapsed 
  # 1.027   0.111   1.139 
system.time(m.sglm <- speedglm(formula=regform, data=testdat, family = binomial()))
m.sglm$coef
  #  user  system elapsed 
  # 0.211   0.055   0.265 
design_mat <- as.matrix(cbind(Intercept=1,testdat[,predvarnms]))
system.time(m.glm.fit <- glm.fit(x=design_mat, y=testdat[,"A"], family = binomial()))
m.glm.fit$coef
  #  user  system elapsed 
  # 0.295   0.081   0.377 
system.time(m.sglm.fit <- speedglm.wfit(y=testdat[,"A"], X=design_mat, family = binomial()))
  #  user  system elapsed 
  # 0.107   0.034   0.151 
m.sglm.fit$coef

fun.glm <- function() {
  testdat <- gen_testdat(predvarnms)
  head(testdat)
  m.glm <- glm(formula = regform, data = testdat, family = binomial())
  m.glm$coef
}
fun.glm.fit <- function() {
  testdat <- gen_testdat(predvarnms)
  m.sglm <- speedglm(formula = regform, data = testdat, family = binomial())
  m.sglm$coef
}
fun.speedglm <- function() {
  testdat <- gen_testdat(predvarnms)
  design_mat <- as.matrix(cbind(Intercept = 1, testdat[,predvarnms]))
  m.glm.fit <- glm.fit(x=design_mat, y=testdat[,"A"], family = binomial())
  m.glm.fit$coef
}
fun.speedglm.wfit <- function() {
  testdat <- gen_testdat(predvarnms)
  design_mat <- as.matrix(cbind(Intercept = 1, testdat[,predvarnms]))
  m.sglm.fit <- speedglm.wfit(y=testdat[,"A"], X=design_mat, family = binomial())
  m.sglm.fit$coef
}
library(microbenchmark)
nsamp <- 10000
speed <- microbenchmark(times=50,
  fun.glm(),
  fun.glm.fit(),
  fun.speedglm(),
  fun.speedglm.wfit()
)
speed
(speed <- mb_summary(speed))
# Unit: milliseconds
#                 expr      min        lq      mean    median        uq      max neval cld
#            fun.glm() 91.64860 102.14260 108.79396 108.28316 114.58223 132.1035    50   c
#        fun.glm.fit() 73.69803  77.99265  83.35009  80.90322  85.59492 105.2169    50  b 
#       fun.speedglm() 77.53957  83.10794  88.76154  85.64683  90.88952 122.4445    50  b 
#  fun.speedglm.wfit() 65.66362  70.34658  76.03511  72.46783  75.37664 188.1003    50 a  
# > (speed <- mb_summary(speed))
#                  name    median
# 1           fun.glm() 108283.16
# 2       fun.glm.fit()  80903.22
# 3      fun.speedglm()  85646.83
# 4 fun.speedglm.wfit()  72467.83
nsamp <- 100000
speed <- microbenchmark(times=10,
  fun.glm(),
  fun.glm.fit(),
  fun.speedglm(),
  fun.speedglm.wfit()
)
mb_summary <- function(x) {
  res <- summary(x, unit="us")
  data.frame(name = res$expr, median = res$median)
}
speed
(speed <- mb_summary(speed))
# Unit: milliseconds
#                 expr       min        lq      mean    median        uq       max neval cld
#            fun.glm() 1164.9370 1226.7263 1268.6743 1264.9778 1306.9116 1364.9334    10   c
#        fun.glm.fit()  817.6741  827.0661  911.7074  888.5717  923.6327 1219.0740    10  b 
#       fun.speedglm()  891.3174  899.4375  952.0213  955.0377  998.3013 1010.7569    10  b 
#  fun.speedglm.wfit()  674.3391  694.1584  761.8763  771.7646  801.1492  878.4608    10 a  
#                  name    median
# 1           fun.glm() 1264977.8
# 2       fun.glm.fit()  888571.7
# 3      fun.speedglm()  955037.7
# 4 fun.speedglm.wfit()  771764.6 






#---------------------------------------------------------------------
### A possible bug in R6. Active binding is called when calling print(obj1) / obj1
# Fixed. Use github version of the package
#---------------------------------------------------------------------
DatNet <- R6Class(classname = "DatNet",
  portable = TRUE,
  class = TRUE,
  public = list(
    names.netVar = character(), 
    dat.netVar = NULL,
    nobs = NA_integer_,        # n of samples in the OBSERVED (original) data
    initialize = function(...) {
      self$nobs <- 5
      self$names.netVar <- c("a", "b", "c")
      self$dat.netVar <- matrix(rnorm(self$nobs*self$netVarcols), nrow = self$nobs, ncol = self$netVarcols)
      invisible(self)
    }
  ),
  active = list(
    netVarcols = function() { length(self$names.netVar) },
    emptydat.netVar = function() {
      self$dat.netVar <- NULL
    }
  ),
  private = list(
    placeholder = list()
  )
)

obj1 <- DatNet$new()
obj1$dat.netVar
obj1

obj1$dat.netVar
obj1






