#' Determine a sensible number of knots
#' 
#' This function computes a reasonable number of inner knots
#' to be used for a basis of unidimensional B-splines
#' @param n An integer vector of sample sizes
#' @return An integer vector with the number of knots for each sample size
#' @references Ruppert, D. (2002). Selecting the number of knots for penalized splines. \emph{Journal of Computational and Graphical Statistics} 11, 735–757.
determine.n.knots <- function(n, cutoff = 35, rate = 0.2) {
#   # Taken from hisemi::n.knots
#   as.integer(trunc(pmin(n, cutoff + pmax(0, n - cutoff)^rate)))
  return(rep(5, length(n)))
}