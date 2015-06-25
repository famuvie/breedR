#' Constructor for a diagonal random effect
#' 
#' Builds a breedR random effect with an index incidence matrix and a diagonal 
#' covariance matrix, with one level for each value of \code{x}.
#' 
#' Uses sparse storage. It does not support nesting (yet). So it is not possible
#' to build random regression coefficients for each level of a grouping factor. 
#' This is in the TODO list.
#' 
#' @param x a numeric covariate or a factor.
#'   
#' @return A random effect with element \code{incidence.matrix} and a diagonal 
#'   \code{structure.matrix}
diagonal <- function (x) {
  
  if (!is.vector(x))
    stopifnot(ncol(as.matrix(x)) == 1)
  x <- as.factor(x)
  
  ## with indMatrix-class, all rows must have exactly one non-zero value
  ## this do not work with diagonal factors with missing values.
  ## Use a more general sparse format that
  if (any(nax <- is.na(x))) {
    inc.mat <- Matrix::sparseMatrix(i = which(!nax),
                                    j = as.integer(x)[!nax],
                                    x = 1,
                                    dims = c(length(x), nlevels(x)))
  } else {
    inc.mat <- as(as.integer(x), 'indMatrix')
  }
  
  colnames(inc.mat) <- levels(x)
  
  eff <- random(incidence = inc.mat,
                covariance = Matrix::Diagonal(nlevels(x)))
  
  class(eff) <- c('diagonal', class(eff))
  
  return(eff)
}
