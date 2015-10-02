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
  
  mc <- match.call()
  arg.list <- as.list(mc)[-1]

  as.Matrix <- function(x) {
    ## If it is not already a Matrix, try conversion to a regular matrix
    ## and then to a Matrix class
    if (!inherits(x, 'Matrix')) {
      x <- as(as.matrix(x), 'Matrix')
    }
    return(x)
  }

  ## Checks ======================================================
  ## arguments
  if (!xor(missing(covariance), missing(precision)))
    stop('Exactly one of covariance and precision matrices
         must be specified in a random effect.')

  for (arg in names(arg.list)) {
    mc[[arg]] <- try(as.Matrix(eval(mc[[arg]], parent.frame())))
    if (inherits(mc[[arg]], 'try-error'))
      stop(paste(arg, 'is expected to be a matrix-like object in a random effects model.'))
  }
  
  ## Define structure either with covariance or precision
  type      <- grep("covariance|precision",
                    names(arg.list),
                    value = TRUE)
  structure <- eval(mc[[type]], parent.frame())
  
  ## conformance
  if (ncol(incidence) != nrow(structure))
    stop('Non conforming incidence and structure matrices.')

  ## Build the random effect, and further specify the generic class
  eff <- breedr_effect(incidence = mc[['incidence']])
  ans <- structure(c(eff,
                     list(structure.matrix = structure,
                          structure.type   = type)),
                   class = c('random', class(eff)))
  
  return(ans)
}
