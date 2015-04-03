#' Get the Pedigree from an object
#' 
#' Returns an object from the formal class \code{pedigree}.
#' @param x object to extract pedigree from
#' @param ... Arguments to be passed to methods.
#' @references \code{\link[pedigreemm]{pedigree-class}} from package
#'   \code{pedigreemm}
#' @export
get_pedigree <- function(x, ...) UseMethod('get_pedigree')

#' Extract the number of traits
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @family metagene
#' @export
get_ntraits <- function(x, ...) UseMethod('get_ntraits')

#' Number of generations
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @family metagene
#' @export
ngenerations <- function(x, ...) UseMethod('ngenerations')

#' Number of individuals
#' @param x a metagene object.
#' @param ... Arguments to be passed to methods.
#' @family metagene
#' @export
nindividuals <- function(x, ...) UseMethod('nindividuals')

#' Simulate a spatial structure
#' 
#' Takes a metagene simulated dataset and distributes its individuals randomly 
#' into a more or less square spatial region. Furthermore, it takes part of the 
#' phenotypic noise and puts some spatial structure into it.
#' 
#' Founders are not put into place, as they don't have phenotypic values. 
#' Therefore, they are returned without coordinates nor spatial values.
#' 
#' The variance of the spatial field is given as a proportion of the variance of
#' the random noise that was added to the Breeding Values to produce the 
#' phenotypes. The phenotypes are afterwards recalculated to preserve the 
#' heritability.
#' 
#' The spatial unit is the distance between consecutive trees.
#' @param meta A metagene object
#' @param variance A number between 0 and 1. The variance of the spatial field 
#'   as a proportion of the non-inheritable pehnotype variance. See Details.
#' @param range A number between 0 and 1. The range of the spatial field, as a 
#'   proportion of the region size.
#' @param ... Arguments to be passed to methods.
#' @return Another metagene object with spatial structure given through an 
#'   additional column \code{sp_X} with the spatially structured component of 
#'   the phenotype, and a 'spatial' list element with coordinates and simulation
#'   information
#' @export
sim.spatial <- function(meta, variance, range, ...) UseMethod('sim.spatial')
