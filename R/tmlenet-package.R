#' Targeted Maximum Likelihood Estimation for Network Data
#'
#' The \pkg{tmlenet} R package implements the Targeted Maximum Likelihood Estimation (TMLE) of causal effects 
#'  under single time point stochastic interventions in network data. 
#'  The package also implements the Horvitz-Thompson estimator for networks (IPTW) and the parametric g-computation formula estimator. 
#'  The inference for the TMLE is based on its asymptotic normality and the efficient influence curve for dependent data. 
#'  The inference for IPTW is based on its corresponding influence curve for dependent data.
#'
#' The input data structure consists of rows of unit-specific observations, with each row \code{i} represented by random variables
#'  (\code{F.i},\code{W.i},\code{A.i},\code{Y.i}), where \code{F.i} is a vector of "\emph{friend IDs}" of unit \code{i} 
#'  (also referred to as \code{i}'s "\emph{network}"), \code{W.i} is a vector of \code{i}'s baseline covariates, \code{A.i} is \code{i}'s exposure
#'  (either binary, categorical or continuous) and \code{Y.i} is \code{i}'s binary outcome.
#'  Each exposure \code{A.i} depends on (possibly multivariate) baseline summary measure(s) \code{sW.i}, where \code{sW.i}
#'  can be any user-specified function of \code{i}'s baseline covariates \code{W.i} and the baseline covariates of \code{i}'s friends in
#'  set \code{F.i} (all \code{W.j} such that \code{j} is in \code{F.i}).
#'  Similarly, each outcome \code{Y.i} depends on \code{sW.i} and (possibly multivariate) summary measure(s) \code{sA.i},
#'  where \code{sA.i} can be any user-specified function of \code{i}'s baseline covariates and exposure (\code{W.i},\code{A.i}) and the
#'  baseline covariates and exposures of \code{i}'s friends (all \code{W.j},\code{A.j} such that \code{j} is in \code{i}'s friend set \code{F.i}). 
#'  
#' The summary measures (\code{sW.i},\code{sA.i}) are defined simultaneously for all \code{i} with functions 
#'  \code{\link{def.sW}} and \code{\link{def.sA}}.
#'  It is assumed that (\code{sW.i},\code{sA.i}) have the same dimensionality across \code{i}. The function \code{\link{eval.summaries}} 
#'  can be used for evaluating these summary measures. All estimation is performed by calling the \code{\link{tmlenet}} function. 
#'  The vector of friends \code{F.i} can be specified either as a single column in the input data (where each \code{F.i} is a string
#'  of friend IDs or friend row numbers delimited by character value in \code{sep}) or as a separate input matrix of network IDs
#'  (where each row is a vector of friend IDs or friend row numbers).
#'  Specifying the network as a matrix generally results in significant improvements to run time.

# @section Documentation:
# \itemize{
# \item To see the package vignette use: \code{vignette("tmlenet_vignette", package="tmlenet")} 
# \item To see all available package documentation use: \code{help(package = 'tmlenet')}
# }
#'
#' @section Routines:
#' The following routines will be generally invoked, in the same order as presented below.
#' \describe{
#' 
#' \item{\code{\link{def.sW}}}{Define the (multivariate) baseline-covariate-based network summary measure functions 
#'    that will be later applied to the input data to obtain a matrix of summary measures \code{sW}. 
#'    Each row \code{i} in \code{sW} will represent a vector of summary measures for observation \code{i}. 
#'    The components (columns) of \code{sW} for unit \code{i} can be defined by any R function that takes as input
#'    \code{i}'s baseline covariates and the baseline covariates of \code{i}'s friends.
#'    This is the first part of the two part specification of the structural equation model for the outcome \code{Y}.
#'    The function \code{def.sW} can take any number of arguments, each argument can specify either a single summary measure
#'    or a multivariate summary measure.
#'    Each argument to \code{def.sW} must be named and has to be a valid R expression or a string containing an R expression.
#'    The expressions can contain double bracket indexing to reference variable values of friends. For example,
#'    \code{Var[[1]]} produces a summary measure that is a vector of length \code{nrow(data)} that contains the values of \code{Var} variable of 
#'    each unit's first friend. One could also use vectors for indexing friends, using a special constant \code{Kmax} to represent the largest number of friends for any unit.
#'    For example, \code{Var[[1:Kmax]]} produces a summary measure that is a matrix of dimension with nrow(data) rows and Kmax columns, 
#'    where the first column will be equal to \code{Var[[1]]}, and so on, up to the last column being equal to \code{Var[[Kmax]]}. Note that \code{Kmax} can vary depending 
#'    on a particular dataset and is determined by the input network structure. That is \code{Kmax} is the largest node degree observed among all observations. 
#'    By default, \code{NA} is used for the missing network summary measure values, however \code{NA}'s can be automatically replaced with 0's 
#'    by passing \code{replaceNAw0 = TRUE} as an argument to \code{def.sW} function.}
#' \item{\code{\link{def.sA}}}{Defines treatment summary measures \code{sA} that can be functions of each unit's exposure & baseline covariates, 
#'    as well the exposures and baseline covariates of unit's friends. 
#'    This is the second part of the two part specification of the structural equation model for the outcome \code{Y}. 
#'    The syntax is identical to \code{def.sW} function, except that \code{def.sA} can consists of functions of baseline covariates
#'    as well as the exposure \code{Anode}.}
#' \item{\code{\link{eval.summaries}}}{A convenience function that may be used for checking and validating the user-defined summary measures. 
#'    Takes the input dataset and evaluates the summary measures based on objects previously defined with function calls \code{def.sW} and \code{def.sA}.
#'    Note that this function is automatically executed by \code{tmlenet} and does not need to be called by the user prior to calling \code{tmlenet}.}
#' \item{\code{\link{tmlenet}}}{Performs estimation of the causal effect of interest using the observed input \code{data}, 
#'    the intervention of interest, the network information and the previously defined summary measures \code{sW}, \code{sA}.}
#' }
#' 
#' @section Datasets:
#' To learn more about the type of data input required by \code{\link{tmlenet}}, see the following example datasets:
#' \itemize{
#'   \item \code{\link{df_netKmax2}}
#'   \item \code{\link{df_netKmax6}}
#'   \item \code{\link{NetInd_mat_Kmax6}}
#' }
#'
#' @section Updates:
#' Check for updates and report bugs at \url{http://github.com/osofr/tmlenet}.
#'
#' @docType package
#' @name tmlenet-package
#'
NULL

#' An example of a row-dependent dataset with known network of at most 2 friends.
#'
#' Simulated dataset containing measured i.i.d. baseline covariate (\code{W1}), dependent binary exposure (\code{A})
#'  and binary binary outcome (\code{Y}), along with a known network of friends encoded by strings on space separated 
#'  friend IDs in \code{Net_str}.
#'  The 1,000 baseline covariates \code{W1} were sampled as i.i.d., 
#'  while the exposure value of \code{A} for each observation \code{i} was sampled 
#'  conditionally on the value of \code{i}'s baseline covariate \code{W1[i]},
#'  as well as the baseline covariate values of \code{i}'s friends in \code{Net_str}. 
#'  Similarly, the binary outcome \code{Y} for each observation was generated conditionally on \code{i}'s 
#'  exposure and baseline covariates values in (\code{W1[i]},\code{A[i]}), 
#'  as well as the values of exposures and baseline covariates of \code{i}'s friends in \code{Net_str}.
#'  Individual variables are described below.
#'
#' @format A data frame with 1,000 dependent observations (rows) and 6 variables:
#' \describe{
#'   \item{IDs}{unique observation identifier}
#'   \item{Y}{binary outcome that depends on unit's baseline covariate value and exposure in \code{W1}, \code{A}, as well as the 
#'      baseline covariate values and exposures \code{W1}, \code{A} of observations in the friend network \code{Net_str}}
#'   \item{nFriends}{number of friends for each observation (row), range 0-2}
#'   \item{W1}{binary baseline covariate (independent)}
#'   \item{A}{binary exposure status that depends on unit's baseline covariate value in \code{W1}, as well as the 
#'      baseline covariate values \code{W1} of observations in the friend network \code{Net_str}}
#'   \item{Net_str}{each observation is a string of space separated friend IDs (this can be either observation IDs or
#'      just space separated friend row numbers)}
#' }
#' @docType data
#' @keywords datasets
#' @name df_netKmax2
#' @usage data(df_netKmax2)
NULL

#' An example of a row-dependent dataset with known network of at most 6 friends.
#'
#' Simulated dataset containing 3 measured i.i.d. baseline covariates (\code{W1}, \code{W2}, \code{W3}), dependent binary exposure (\code{A})
#'  and dependent binary binary outcome (\code{Y}), along with a known network of friends encoded by strings on space separated
#'  friend IDs in \code{Net_str}.
#'  The baseline covariates (\code{W1},\code{W2},\code{W3}) were sampled as i.i.d.,
#'  while the exposure value of \code{A} for each observation \code{i} was sampled
#'  conditionally on the values of \code{i}'s baseline covariates (\code{W1[i]} \code{W2[i]}, \code{W3[i]}),
#'  as well as the baseline covariate values of \code{i}'s friends in \code{Net_str}.
#'  Similarly, the binary outcome \code{Y} for each observation was generated conditionally on \code{i}'s
#'  exposure and baseline covariates values in (\code{W1[i]},\code{W2[i]},\code{W3[i]},\code{A[i]}),
#'  as well as the values of exposures and baseline covariates of \code{i}'s friends in \code{Net_str}.
#'  Individual variables are described below.
#'
#' @format A data frame with 1,000 dependent observations (rows) and 6 variables:
#' \describe{
#'   \item{IDs}{unique observation identifier}
#'   \item{W1}{categorical baseline covariate (independent), range 0-5}
#'   \item{W2}{binary baseline covariate (independent)}
#'   \item{W3}{binary baseline covariate (independent)}
#'   \item{A}{binary exposure that depends on unit's baseline covariate values, as well as the
#'      baseline covariate values of observations in the friend network \code{Net_str}}
#'   \item{Y}{binary outcome that depends on unit's baseline covariate value and exposure, as well as the
#'      baseline covariate values and exposures of observations in the friend network \code{Net_str}}
#'   \item{nFriends}{number of friends for each observation (row), range 0-6}
#'   \item{Net_str}{a vector of strings, where for each observation its a string of space separated friend IDs (this can
#'      be either observation IDs or just space separated friend row numbers)}
#' }
#' @docType data
#' @keywords datasets
#' @name df_netKmax6
#' @usage data(df_netKmax6)
NULL

#' An example of a network ID matrix
#'
#' This matrix demonstrates an internal data structure of how \code{tmlenet} stores the input network.
#' The network in this matrix is derived from the column \code{"Net_str"} of the dataset \code{df_netKmax6}.
#' Both, this matrix and the vector of strings in column \code{"Net_str"}, represent the same network and both can be used 
#' for specifying the input network to \code{tmlenet} function.
#' The network matrix may be specified by using the argument \code{NETIDmat} of the \code{tmlenet} function.
#' Inputing the network via this type of a matrix may lead to significant reduction in total run time,
#' since any network specified as a vector of strings, such as in column "Net_str",
#' will be first converted to this type of matrix representation.
#' 
#' See below and Example 3 in \code{tmlenet} help file for examples constructing this matrix from the initial network input in 
#' column \code{"Net_str"} of \code{df_netKmax6}.
#' 
#' This matrix consists of \code{1000} rows and \code{6} columns. Each row \code{i} encodes a vector of network IDs of observation \code{i} in 
#' \code{df_netKmax6}, i.e., 
#' \code{NetInd_mat_Kmax6[,i]} contains a vector of observation row numbers in \code{df_netKmax6} that are presumed "connected to" (or "friends of") 
#' observation \code{i}. Each observation can have at most 6 friends and if an observation \code{i} has fewer than 6 friends the remainder row 
#' of \code{NetInd_mat_Kmax6[,i]} is filled with \code{NA}s.
#' @docType data
#' @keywords datasets
#' @name NetInd_mat_Kmax6
#' @usage data(NetInd_mat_Kmax6)
#' @examples
#' 
#'data(df_netKmax6)
#'Net_str <- df_netKmax6[, "Net_str"]
#'IDs_str <- df_netKmax6[, "IDs"]
#'net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = ncol(df_netKmax6))
#'net_ind_obj$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = ' ')
#'NetInd_mat <- net_ind_obj$NetInd
NULL







