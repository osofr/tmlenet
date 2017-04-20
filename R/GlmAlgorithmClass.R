#' logisfitR6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
logisfitR6 <- R6Class("logisfitR6",
  public =
    list(
         lmclass = NULL,
         fitfunName = NULL,
         initialize = function() {
           self$is_valid
         },

         fit = function(Xmat, Yvals) {
           stop('Override this function in a subclass')
         },
 
         predict = function() {
           stop('Override this function in a subclass')
         },

         update = function() {
           stop('Override this function in a subclass')
         },
 
         process = function(datsum_obj) {
           if (gvars$verbose) print(paste("calling glm.generic for", self$fitfunname))
           Xmat <- datsum_obj$getXmat
           Y_vals <- datsum_obj$getY
 
           # Xmat has 0 rows: return NA's and avoid throwing exception:
           if (nrow(Xmat) == 0L) {
             m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
           } else {
             m.fit <- self$fit(Xmat, Y_vals)
           }
           fit <- list(coef = m.fit$coef, linkfun = "logit_linkinv", fitfunname = self$fitfunname)
           if (gvars$verbose) print(fit$coef)
           class(fit) <- c(class(fit), c(self$lmclass))
           return(fit)
         }
         ),
  active =
    list(
         is_valid = function() {
           errors = character()
           if(is.null(self$fitfunname)){
             errors <- c(errors, 'Please define a fitfunname')
           }
           if(is.null(self$lmclass)){
             errors <- c(errors, 'Please define a lmclass')
           }

           # Defines fit
           # Defines predict
           if(length(errors)!= 0) stop(errors)
           TRUE
         }
        ),
  private =
    list(
    )
)

#' glmR6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
glmR6 <- R6Class("glmR6",
  inherit = logisfitR6,
  public =
    list(
        fitfunname='glm',
        lmclass='glmR6',

        initialize = function() {
        },

        fit = function(Xmat, Yvals) {
          ctrl <- glm.control(trace = FALSE)
          # ctrl <- glm.control(trace = FALSE, maxit = 1000)
          SuppressGivenWarnings({
            return(stats::glm.fit(x = Xmat, y = Yvals, family = binomial() , control = ctrl))
          }, GetWarningsToSuppress())
        },

        predict = function() {
        }

        ),
  active =
    list(
        ),
  private =
    list(
    )
)

#' speedglmR6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
speedglmR6 <- R6Class("speedglmR6",
  inherit = logisfitR6,
  public =
    list(
        fitfunname='speedglm',
        lmclass='speedglmS3',

        initialize = function() {
        },

        fit = function(Xmat, Yvals) {
          # , maxit=1000
          m.fit <- try(speedglm::speedglm.wfit(X = Xmat, y = Yvals, family = binomial(), trace = FALSE, method='Cholesky'), silent = TRUE)
          if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
            message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", m.fit)
            return(private$fallback_function(Xmat, Yvals))
          }
          return(m.fit)
        },

        predict = function() {
          
        }
        ),
  active =
    list(
        ),
  private =
    list(
         fallback_function = glmR6$new()$fit
    )
)
