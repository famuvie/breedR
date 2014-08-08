## Random regression models

#' Evaluate Legendre polynomials
#' 
#' Evaluate the Legendre polynomials up to certain order on a given set
#' of values. Normalization is done within the function.
#' 
#' @param x numeric. Vector of values for evaluation after normalization into [-1, 1]
#' @param ord integer. Maximum order of the polynomials to evaluate.
#' 
#' @export
legendre_values <- function (x, ord) {
  # Normalize x
  rg <- range(x)
  w <- 2*(x - rg[1])/diff(rg) - 1
  
  # Legendre polynomials list up to order ord
  p.list <- legendre.polynomials(ord, normalized = TRUE)
  
  # Polynomial values for each order
  v.list <- polynomial.values(p.list, w)
  
  # Return matrix-wise
  return(do.call(cbind, v.list))
}
