#' Constructor for a random effect
#' 
#' The random-class is virtual. No object should be directly created with this 
#' constructor. This constructor is to be called from within subclasses like 
#' generic, diagonal, spatial or additive_genetic.
#' 
#' This constructor performs the arguments and conformance checks. But the
#' implementation details (i.e., storage format and handling) is left for the
#' subclasses.
#' 
#' @param incidence matrix-like object
#' @param covariance matrix-like object
#' @param precision matrix-like object
#'   
#' @return A list with elements \code{incidence.matrix}, \code{structure.matrix}
#'   and \code{structure.type}, which is a string indicating either 
#'   \code{covariance} or \code{precision}.
random <- function (incidence, covariance, precision) {
  
  ## Checks ======================================================
  ## arguments
  if (!xor(missing(covariance), missing(precision)))
    stop('Exactly one of covariance and precision matrices
         must be specified in a random effect.')

  ## Define structure either with covariance or precision
  mc <- match.call()
  arg.list <- as.list(mc)[-1]
  type      <- grep("covariance|precision",
                    names(arg.list),
                    value = TRUE)
  structure <- eval(arg.list[[type]])
  
  ## conformance
  if (ncol(incidence) != nrow(structure))
    stop('Non conforming incidence and structure matrices.')

  ## Build the random effect, and further specify the generic class
  eff <- breedr_effect(incidence = incidence)
  ans <- structure(c(eff,
                     list(structure.matrix = structure,
                          structure.type   = type)),
                   class = c('random', class(eff)))
  
  return(ans)
}


#' @importFrom stats vcov
vcov.random <- function(object, ...) {
  ans <- object$structure.matrix
  attr(ans, 'inverse') <- object$structure.type == 'precision'
  ans
}
