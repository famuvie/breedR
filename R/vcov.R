
#' @importFrom stats vcov
#' @export
vcov.random <- function(object, ...) {
  ans <- object$structure.matrix
  attr(ans, 'inverse') <- object$structure.type == 'precision'
  ans
}
