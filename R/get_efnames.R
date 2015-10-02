
#' Get names of effects
#' 
#' Give the names of the (components) of the effects
#' in a breedr_modelframe. Internal function. Not exported.
#' 
#' @param effects A breedr_modelframe.
get_efnames <- function(effects) {
  
  subnames <- function(idx) {
    sizes <- vapply(effects, effect_size, 1)
    if (sizes[idx] <= 1)
      return(names(effects)[idx]) else
        return(names(effects[[idx]]$effects))
  }
  ans <- unlist(lapply(seq_along(effects), subnames))
  return(ans)
}
