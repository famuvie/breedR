#' Constructor for a generic effect
#' 
#' The breedr_effect-class is virtual. No object should be directly created with
#' this constructor. This constructor is to be called from within non-virtual 
#' subclasses like generic, diagonal, spatial or additive_genetic.
#' 
#' This constructor performs the arguments checks. But the implementation
#' details (i.e., storage format and handling) is left for the subclasses.
#' 
#' @param incidence matrix-like object
#'   
#' @return A list with a single element \code{incidence.matrix}.
breedr_effect <- function(incidence) {
  
  ## Checks ======================================================
  ## arguments
  if (missing(incidence) )
    stop('Incidence matrix must be specified in a random-effects model.')

  ## Build the random effect, and further specify the generic class
  ans <- list(incidence.matrix = incidence)
  class(ans) <- c('breedr_effect')
  
  return(ans)
}

#' Constructor for a group of effects
#' 
#' Builds an \code{effect_group} from a list of \code{breer_effect} elements.
#' 
#' Temporarily, this takes the \code{cov.ini} argument and includes it in the 
#' object. In the future, the initial covariance matrix will be a matter of the
#' inference engine, not inherent to the model.
#' 
#' @param x list of breedr_effect elements
#' @param cov.ini initial covariance matrix for the estimation algorithm
#' 
#' @return A list of \code{breedr_effect} elements.
effect_group <- function(x, cov.ini) {
  
  ## Checks ==========================================
  ## x is a list and cov.ini a SPD matrix
  stopifnot(is.list(x))
  cov.ini <- as.matrix(cov.ini)
  stopifnot(isSymmetric.matrix(cov.ini))
  ev <- eigen(cov.ini, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(ev > 0))
  
  ## all elements are breer_effects
  if (!all(sapply(x, inherits, 'breedr_effect')))
    stop('All of the effects must be of class breedr_effect.')
  
  ## cov.ini is square and of size equal to number of effects
  nx <- length(x)
  if (!all(dim(cov.ini) == nx))
    stop('Dimension of the initial covariance matrix do not conform with
         number of effects in the group.')
  
  ans <- structure(list(effects = x, 
                        cov.ini = cov.ini),
                   class = 'effect_group')
  return(ans)
}

#' Size of a group of effects
#' 
#' @param x object of class \code{effect_group}
group_size <- function(x) {
  stopifnot(inherits(x, 'effect_group'))
  length(x$effects)
}
