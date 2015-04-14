#' Build a generic model
#' 
#' Check conformity of arguments and return a \code{generic} object.
#' 
#' A generic random effect stores the incidence and structure matrices in Matrix
#' form, which tries to take advantage of sparsity, if it exists.
#' 
#' @param incidence matrix-like object
#' @param covariance matrix-like object
#' @param precision matrix-like object
#'   
#' @return A list with elements \code{incidence.matrix}, \code{structure.matrix}
#'   and \code{structure.type}, which is a string indicating either 
#'   \code{covariance} or \code{precision}.
generic <- function(incidence, covariance, precision) {
  
  ## Try converting to Matrix class (either sparse or dense)
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
  for (arg in names(arg.list)) {
    mc[[arg]] <- try(as.Matrix(eval(mc[[arg]], parent.frame())))
    if (inherits(mc[[arg]], 'try-error'))
      stop(paste(arg, 'is expected to be a matrix-like object in a generic model.'))
  }

  ## Build the random effect, and further specify the generic class
  random.call <- mc
  random.call[[1]] <- as.symbol('random')
  ans <- eval(random.call)
  class(ans) <- c('generic', class(ans))
  
  return(ans)
}

