#' @describeIn get_param Get the param from a remlf90 object
get_param.remlf90 <- function(x) {
  
  return(get_param.breedr_modelframe(x$effects))
  
}

#' @describeIn get_param Get the param from a breedr_modelframe object.
#'   Internal function.
get_param.breedr_modelframe <- function(x) {
  
  ## Identify effect_groups
  eg.idx <- which(vapply(x, inherits, TRUE, 'effect_group'))
  
  ## If no param, return NULL
  if (length(eg.idx) == 0) return(NULL)
  
  param.list <- lapply(x[eg.idx], get_param.effect_group)
  
  return(param.list)
}


#' @describeIn get_param Get the param from a effect_group object.
#'   Internal function.
get_param.effect_group <- function(x) {
  
  ## Identify spatial elements
  sp.idx <- which(vapply(x$effects, inherits, TRUE, 'spatial'))
  
  ## If no spatial elements, return NULL
  if (length(sp.idx) == 0) return(NULL)
  
  ## Otherwise, there must be at least one param
  param.list <- lapply(x$effects[sp.idx], get_param.spatial)
  stopifnot(length(param.list) > 0)

  ## If more than one param in the group, all must be identical
  if ((np <- length(param.list)) > 1) {
    for (i in seq_len(np - 1)) {
      stopifnot(identical(param.list[[1]], param.list[[i+1]]))
    }
  }
  
  return(param.list[[1]])
}

#' @describeIn get_param Get the param from a \code{spatial} object
#' @family spatial
#' @export
get_param.spatial <- function(x) {
  return(x$param)
}
