#' @describeIn get_pedigree Get the pedigree from a \code{metagene} object
#' @family metagene
#' @export
get_pedigree.metagene <- function(x, ...) {
  return(with(x$Data, pedigreemm::pedigree(sire=dad, dam=mum, label=self)))
}


#' @describeIn get_pedigree Get the pedigree from a remlf90 object
#' @export
get_pedigree.remlf90 <- function(x, ...) {
  
  return(get_pedigree.breedr_modelframe(x$effects))
  
  #   ped <- x$effects$genetic$ped
  #   if( !is.null(ped) ) {
  #     map <- attr(ped, 'map')
  #     ped <- with(ped,
  #                 pedigreemm::pedigree(sire=sire, dam=dam, label=self))
  #     attr(ped, 'map') <- map
  #   }
  #   return(ped)
}

#' @describeIn get_pedigree Get the pedigree from a breedr_modelframe object.
#'   Internal function.
get_pedigree.breedr_modelframe <- function(x, ...) {
  
  ## Identify effect_groups
  eg.idx <- which(vapply(x, inherits, TRUE, 'effect_group'))
  
  ## If no effect_group (i.e. no pedigree), return NULL
  if (length(eg.idx) == 0) return(NULL)
  
  ped.list <- lapply(x[eg.idx], get_pedigree.effect_group)
  ped.list <- ped.list[!vapply(ped.list, is.null, TRUE)]
  if (length(ped.list) == 1) ped.list <- ped.list[[1]]

  ## If no pedigree, return NULL
  if (length(ped.list) == 0) return(NULL)
  
  return(ped.list)
}


#' @describeIn get_pedigree Get the pedigree from a effect_group object.
#'   Internal function.
get_pedigree.effect_group <- function(x, ...) {
  
  ## Identify genetic elements
  ge.idx <- which(vapply(x$effects, inherits, TRUE, 'genetic'))
  
  ## If no genetic elements, return NULL
  if (length(ge.idx) == 0) return(NULL)
  
  ## Otherwise, there must be at least one pedigree
  ped.list <- lapply(x$effects[ge.idx], get_pedigree.genetic)
  stopifnot(length(ped.list) > 0)

  ## If more than one pedigree in the group, all must be identical
  if ((np <- length(ped.list)) > 1) {
    for (i in seq_len(np - 1)) {
      stopifnot(identical(ped.list[[1]], ped.list[[i+1]]))
    }
  }
  
  return(ped.list[[1]])
}

#' @describeIn get_pedigree Get the pedigree from a \code{genetic} object
#' @family genetic
#' @export
get_pedigree.genetic <- function(x, ...) {
  return(x$pedigree)
}
