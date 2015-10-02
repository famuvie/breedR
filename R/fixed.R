#' Constructor for a fixed effect
#' 
#' @param x a numeric covariate or a factor.
#'   
#' @return A breedr_effect with element \code{incidence.matrix}.
fixed <- function (x) {
  
  eff <- breedr_effect(incidence = x)
  class(eff) <- c('fixed', class(eff))
  
  return(eff)
}
