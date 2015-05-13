## Functions for checking model components
## Internal - not exported

## TODO: 
##  - check_spatial
##  - check_generic
##  - make tests for them

check_genetic <- function(model = c('add_animal', 'competition'),
                          pedigree,
                          id,
                          coordinates,
                          competition_decay = 1,
                          pec = list(present = FALSE),
                          autofill = TRUE,
                          var.ini) {
  mc <- match.call()
  
  for (arg in c('model', 'pedigree', 'id', 'var.ini')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required in the genetic component.'))
  }
  
  if (!inherits(pedigree, c('pedigree','data.frame'))){
    stop(paste('Argument', pedigree, 'must be an object of type pedigree or data.frame'))
  }
  
  mc$model <- match.arg(model)
  
  if(mc$model == 'add_animal'){
    if (var.ini <= 0 || !is.numeric(var.ini) )
     stop (paste('var.ini must be a positive number'))
  }
  
  if (mc$model == 'competition') {
    for (arg in c('coordinates')) {
      if (eval(call('missing', as.name(arg))))
        stop(paste('Argument', arg, 'required in the genetic component.'))
    }
    
    mc$coordinates <- try(as.data.frame(coordinates))
    if (inherits(mc$coordinates, 'try-error'))
      stop('Argument coordinates in the genetic component not coercible to a data.frame')
    if (!all(vapply(coordinates, is.numeric, TRUE)))
      stop('Argument coordinates in the genetic component not numeric')
    if (ncol(coordinates) != 2)
      stop('Only two dimensions admitted for coordinates in the genetic component.')
    if (!all(vapply(pec, is.logical,TRUE)) || !is.list(pec) || !'names' %in% names(attributes(pec)))
     stop('pec must be a named list with logical elements')
    if (nrow(var.ini)!=ncol(var.ini))
      stop('var.ini must be a square matrix')
    if(!all(var.ini == t(var.ini)))
      stop('var.ini must be a SPD matrix')
    if (!all(eigen(var.ini)$values >0))
     stop('var.ini must be a SPD matrix')
  }
  
  stopifnot(is.numeric(competition_decay))
  stopifnot(competition_decay > 0)
  
  return(as.list(mc[-1]))
}



check_spatial <- function(model = c('splines', 'AR'),
                          coordinates,
                          n.knots = c(6,6) ,
                          rho = c(0.5,0.5),
                          var.ini) {
  mc <- match.call()
    
  for (arg in c('model', 'coordinates', 'var.ini')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required in the spatial component.'))
  }
  
  mc$coordinates <- try(as.data.frame(coordinates))
  
  if (inherits(mc$coordinates, 'try-error'))
    stop('Argument coordinates in the spatial component not coercible to a data.frame')
  if (!all(vapply(coordinates, is.numeric, TRUE)))
    stop('Argument coordinates in the spatial component is not numeric')
  if (ncol(coordinates) != 2)
    stop('Only two dimensions admitted for coordinates in the spatial component.')
  
  mc$model <- match.arg(model)
  
  if(mc$model == 'splines'){
    if (!is.vector(n.knots) || length(n.knots) !=2 || !all(n.knots%%1==0))
      stop(paste('n.knots must be a vector of two integers'))
  }
  
  if (mc$model == 'AR'){ # To do : Check if rho is a matrix or data frame.
    if(is.vector(rho)){
      if (length(rho)!=2)
        stop('rho must contain exactly two components')
      if (!all(vapply(rho, is.numeric, TRUE)))
        stop('Argument rho in the spatial component is not numeric')
      if (any(abs(rho)>=1))
        stop('rho must contain two numbers strictly between -1 and 1')
    }
  }
  
  if (var.ini <= 0 || is.numeric(var.ini) == FALSE )
    stop (paste('var.ini must be a positive number'))
  
  return(as.list(mc[-1]))
}



check_generic <- function(){
  
}

