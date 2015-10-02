#' Build a generic model
#' 
#' Check conformity of arguments and return a \code{generic} object.
#' 
#' A generic random effect stores the incidence and structure matrices in Matrix
#' form, which tries to take advantage of sparsity, if it exists.
#' 
#' @inheritParams random
#' @return A list with elements \code{incidence.matrix}, \code{structure.matrix}
#'   and \code{structure.type}, which is a string indicating either 
#'   \code{covariance} or \code{precision}.
generic <- function(incidence, covariance, precision) {
  
  ## Try converting to Matrix class (either sparse or dense)
  mc <- match.call()

  ## Build the random effect, and further specify the generic class
  random.call <- mc
  random.call[[1]] <- as.symbol('random')
  ans <- eval(random.call, parent.frame())
  class(ans) <- c('generic', class(ans))
  
  return(ans)
}

