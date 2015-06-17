#' @describeIn get_pedigree Get the pedigree from a \code{metagene} object
#' @family metagene
#' @export
get_pedigree.metagene <- function(x, ...) {
  return(with(x$Data, pedigreemm::pedigree(sire=dad, dam=mum, label=self)))
}


#' @describeIn get_pedigree Get the pedigree from a remlf90 object
#' @export
get_pedigree.remlf90 <- function(x, ...) {
  ped <- x$effects$genetic$ped
  if( !is.null(ped) ) {
    map <- attr(ped, 'map')
    ped <- with(ped,
                pedigreemm::pedigree(sire=sire, dam=dam, label=self))
    attr(ped, 'map') <- map
  }
  return(ped)
}


#' @describeIn get_pedigree Get the pedigree from a \code{genetic} object
#' @family genetic
#' @export
get_pedigree.genetic <- function(x, ...) {
  return(x$pedigree)
}
