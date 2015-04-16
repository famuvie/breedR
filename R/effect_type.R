
#' @describeIn effect_type Type of an \code{effect_group}
#' @export
effect_type.effect_group <- function(x) {
  
  types <- sapply(x$effects, effect_type.breedr_effect)
  utypes <- unique(types)
  stopifnot(length(utypes) == 1)
  return(utypes)
}

#' @describeIn effect_type Type of an \code{breedr_effect}
#' @export
effect_type.breedr_effect <- function(x) {
  ifelse(inherits(x, 'random'), 'random', 'fixed')
}
