## Random regression models

#' Evaluate Legendre polynomials
#' 
#' Evaluate the Legendre polynomials up to certain order on a given set of 
#' values. Normalization is done within the function.
#' 
#' By default,
#' 
#' @param x numeric. Vector of values for evaluation after normalization into 
#'   [-1, 1]
#' @param ord integer. Maximum order of the polynomials to evaluate.
#' @param scale vector of length 2. Scale to use for normalization. See Details.
#'   
#'   By default, the empirical range of \code{x} is linearly scaled into 
#'   [-1, 1]. However, sometimes it is useful to give a predefined range. For
#'   example, when the legendre values are required at several different sets of
#'   points x but in the same consistent scale. Note that \code{range = 0:1}
#'   implies no re-scaling at all.
#' @export
legendre_values <- function (x, ord, scale = range(x)) {
  # Normalize x
  if( diff(scale) > 0 ) {
    w <- 2*(x - scale[1])/diff(scale) - 1
  } else {
    # constant vector or one value only
    w <- 0
  }
  # Legendre polynomials list up to order ord
  p.list <- legendre.polynomials(ord, normalized = TRUE)
  
  # Polynomial values for each order
  v.list <- polynomial.values(p.list, w)
  
  # Return matrix-wise
  return(do.call(cbind, v.list))
}
