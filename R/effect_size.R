
#' @describeIn effect_size Size of an \code{effect_group}
#' @export
effect_size.effect_group <- function(x) {
  
  return(nrow(as.matrix(x$cov.ini)))
}

#' @describeIn effect_size Size of an \code{breedr_effect}
#' @export
effect_size.breedr_effect <- function(x) {
  ifelse(inherits(x, 'random'), 1, 0)
}
