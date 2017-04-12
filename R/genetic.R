#' Build an genetic model
#' 
#' Check conformity of arguments and return a \code{genetic} object.
#' 
#' This is a virtual class. No objects are expected to be created directly.
#' 
#' @param pedigree object of class 'pedigree'
#' @inheritParams random
#' @return A list with elements \code{pedigree}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
genetic <- function(pedigree, incidence, covariance, precision) {
  
  ## checks
  stopifnot(!missing(pedigree))
  if (nrow(as.data.frame(pedigree)) != ncol(incidence))
    stop('The incidence matrix should have as many columns as individuals in the pedigree.')
  
  ## Build the random effect, and further specify the genetic class
  random.call <- mc <- match.call()
  arg.list <- as.list(mc)[-1]
  
  random.call[[1]] <- as.symbol('random')
  ans <- eval(random.call[-match('pedigree', as.list(random.call))],
              parent.frame())

  ans$pedigree <- pedigree
  class(ans) <- c('genetic', class(ans))
  
  return(ans)
}


#' Build an additive_genetic model
#' 
#' Check conformity of arguments and return a \code{additive_genetic} object.
#' 
#' @param pedigree object of class 'pedigree'
#' @inheritParams random
#' @return A list with elements \code{pedigree}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' ped <- pedigreemm::pedigree(sire = c(NA,NA,1, 1,4,5),
#'                             dam  = c(NA,NA,2,NA,3,2),
#'                             label= 1:6)
#' inc <- cbind(0, 0, diag(4))
#' breedR:::additive_genetic(ped, inc)
additive_genetic <- function(pedigree, incidence) {
  
  
  ## Build the genetic effect, and further specify the additive_genetic class
  relationship.matrix <- pedigreemm::getA(pedigree)
  
  ## NOTE:
  ## We could as well use pedigreemm::getAInv(pedigree)
  ## to construct the precision matrix.
  
  ans <- genetic(pedigree, incidence, covariance = relationship.matrix)

  class(ans) <- c('additive_genetic', class(ans))
  
  return(ans)
}


#' Build an additive-genetic animal model
#' 
#' Given a pedigree, and an index vector of observations, build and 
#' \code{additive_genetic_animal} model.
#' 
#' \code{idx} must hold the index of observed individuals in the original 
#' codification. If recoding took place when building the pedigree, this
#' function will convert the codes internally.
#' 
#' @param idx integer vector of observed individuals (in the original
#'   codification)
#' @inheritParams additive_genetic
#' @importFrom methods as
#'   
#' @return A list with elements \code{pedigree}, \code{incidence.matrix}, 
#'   \code{structure.matrix} and \code{structure.type}, which is a string 
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' dat <- data.frame(id = 1:4,
#'                   sire = c(11, 11, 2, 3),
#'                   dam  = c(12, NA, 1, 12))
#' ped <- build_pedigree(1:3, data = dat)
#' breedR:::additive_genetic_animal(ped, dat$id)
additive_genetic_animal <- function(pedigree, idx) {
  
  ## Checks
  stopifnot(is.numeric(idx))
  # Not necessarily: might be multiple observations
  # stopifnot(length(idx) < nrow(as.data.frame(pedigree)))
  
  ## Incidence matrix
  ## It is possible that the pedigree has been recoded/reordered
  ## In that case, we need to recode the data file id codes as well
  if (!is.null(attr(pedigree, 'map')))
    idx <- attr(pedigree, 'map')[idx]
  
  ## All the (recoded) codes in idx must be within pedigree@label
  ## it is enough to check max(idx) since pedigree@label = 1, ..., n
  stopifnot(max(idx) <= length(pedigree@label))
  
  ## The pedigree might potentially have further individuals
  ## to evaluate (either founders, or descendants).
  inc.mat <- as(
    Matrix::sparseMatrix(i = seq_along(idx),
                         j = idx,
                         x = 1,
                         dims = c(length(idx),
                                  nrow(as.data.frame(pedigree)))),
    'indMatrix')
  
  ans <- additive_genetic(pedigree, inc.mat)
  class(ans) <- c('additive_genetic_animal', class(ans))
  return(ans)
}



#' Build an additive-genetic competition model
#' 
#' Return incidence and structure for a \code{additive_genetic_competition}
#' model, given the pedigree, the spatial coordinates and codes of the
#' observations and the competition decay parameter.
#' 
#' \code{id} must hold the codes of observed individuals in the original 
#' codification. If recoding took place when building the pedigree, this 
#' function will handle the codes internally.
#' 
#' @param id integer vector of numeric codes for observed individuals
#' @inheritParams additive_genetic
#' @inheritParams competition
#' @inheritParams build_grid
#'   
#' @return A list with elements \code{pedigree}, \code{incidence.matrix}, 
#'   \code{structure.matrix} and \code{structure.type}, which is a string 
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' dat <- data.frame(id   = 1:5,
#'                   sire = c(11, 11, 2, 3, 2),
#'                   dam  = c(12, NA, 1, 12, 1),
#'                   x    = c(rep(1:2, times = 2), 3),
#'                   y    = c(rep(1:2, each = 2), 3))
#' ped <- build_pedigree(1:3, data = dat)
#' breedR:::additive_genetic_competition(ped, coord = dat[, c('x', 'y')], dat$id, 2)
additive_genetic_competition <- function(pedigree,
                                         coordinates,
                                         id,
                                         decay,
                                         autofill = TRUE) {
  
  ## Checks
  stopifnot(is.numeric(id))
  stopifnot(length(id) == nrow(coordinates))
  # Not necessarily: might be multiple observations per genotype
  # stopifnot( (n <- length(id)) < (p <- nrow(as.data.frame(pedigree))))
  
  n <- length(id)    # n observations
  p <- nrow(as.data.frame(pedigree))   # n genotypes
  
  ## additive_genetic_competition inherits from additive_genetic
  ## and from competition simultaneously.
  ## As S3 classes do not support multiple inheritance, simulate it
  ## manually by constructing separate competition and additive_genetic
  ## objects with dummy covariance and incidence matrices respectively
  ## and then composing manually the random effect
  cov.dummy <- Matrix::Diagonal(n)
  comp.aux <- competition(coordinates = coordinates, 
                          covariance  = cov.dummy,
                          decay       = decay,
                          autofill    = autofill)
  
  inc.dummy <- Matrix::Diagonal(p)
  ag.aux <- additive_genetic(pedigree, inc.dummy)
  
  ## the internal codes of the observed individuals are the indices
  ## of the corresponding levels of the random effect
  if ('map' %in% names(attributes(pedigree))) {
    id.internal  <- attr(pedigree, 'map')[id]
  } else {
    id.internal <- id
  }
  ## the incidence matrix has to be exapanded with 0
  ## in the columns corresponding to unobserved individuals (e.g. founders)
  ## furthermore, the order of the columns must fit
  ## the internal representation
  inc.mat <- Matrix::Matrix(0, n, p)
  inc.mat[, id.internal] <- as.matrix(comp.aux$incidence.matrix)
  
  random.args <- structure(list(inc.mat, ag.aux$structure.matrix),
                           names = c('incidence', ag.aux$structure.type))
  
  ans <- do.call('random', random.args)
  ans$pedigree <- pedigree
  ans$coordinates <- coordinates
  
  ## Include inherited classes from the additive_genetic object 
  ## in the first place
  class(ans) <- c('additive_genetic_competition',
                  setdiff(class(ag.aux), class(ans)),
                  setdiff(class(comp.aux), class(ans)),
                  class(ans))
  
  return(ans)
}

