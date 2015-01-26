#' Get the Pedigree from an object
#' 
#' Returns an object from the formal class \code{pedigree}.
#' @param x object to extract pedigree from
#' @param ... Arguments to be passed to methods.
#' @references \code{\link[pedigreemm]{pedigree-class}} from package
#'   \code{pedigreemm}
#' @family metagene
#' @export
get_pedigree <- function(x, ...) UseMethod('get_pedigree')
