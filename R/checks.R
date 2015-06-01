## Functions for checking model components
## Internal - not exported

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
    if (!is.list(pec))
      stop('pec must be a list')
    if (is.null(names(pec)) || !all(nchar(names(pec))>0))
     stop('pec must be a named list')
    
    if (!setequal(names(pec), c('present', 'var.ini')))
      stop('Unrecognized argument in pec')
    if (!is.logical(pec$present))
      stop('logical value expected in pec$present')
    
    if (nrow(var.ini)!=ncol(var.ini))
      stop('var.ini must be a square matrix')
    if(!all(var.ini == t(var.ini)))
      stop('var.ini must be a SPD matrix')
    if (!all(eigen(var.ini)$values >0))
     stop('var.ini must be a SPD matrix')
    mc$pec <- pec
    mc$competition_decay <- competition_decay
  }
  
  stopifnot(is.numeric(competition_decay))
  stopifnot(competition_decay > 0)
  
  mc$var.ini <- var.ini
  mc$pedigree <- pedigree
  mc$id <- id
  mc$autofill <- autofill
  
  
  return(as.list(mc[-1]))
}



check_spatial <- function(model = c('splines', 'AR'),
                          coordinates,
                          n.knots,
                          rho,
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
    mc$n.knots <- n.knots
  }
  
  if (mc$model == 'AR'){ 
    if (!is.vector(rho))
      stop('rho must be a vector')
    if (length(rho)!=2)
      stop('rho must contain exactly two components')
    if (!all(vapply(rho, is.numeric, TRUE)))
      stop('Argument rho in the spatial component is not numeric')
    if (any(abs(rho)>=1))
      stop('rho must contain two numbers strictly between -1 and 1')
    
    mc$rho <- rho
  }
  
  if (var.ini <= 0 || !is.numeric(var.ini) || length(var.ini) != 1)
    stop (paste('var.ini must be a positive number'))
  
  mc$var.ini <- var.ini
  
  return(as.list(mc[-1]))
}



check_generic <- function(x){
  
  mc <- match.call()
  
  for (arg in c('x')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required in the generic component.'))
  }
  
  if (!is.list(x))
    stop('Argument x in the generic component must be a list')
  if (is.null(names(x)))
    stop('Argument x must be a named list')
  if (!all(nchar(names(x))>0))
    stop('All elements of the argument x must be named')
  if (any(duplicated(names(x))))
    stop('Argument x must be a named list with different names')
  if (!all(sapply(x,is.list)))
    stop('All elements of the argument x must be list elements')
  
  for (arg.idx in seq.int(x)){ 
    result <- try(do.call('valid_generic_element', x[[arg.idx]]), silent =TRUE)
    if (inherits(result, 'try-error'))
      stop(paste(attr(result, 'condition')$message, 'in generic component', names(x)[arg.idx]))
  }
  
  mc$x <- x
  
  return(as.list(mc[-1]))
}


valid_generic_element <- function(incidence, covariance, precision, var.ini){
  
  for (arg in c('incidence', 'var.ini')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required'))
  }
  if (!xor(missing(covariance), missing(precision)))
    stop(paste('Exactly one argument between covariance and precision must be specified'))
  
  if (missing(covariance)) {
    structure <- precision
    str.name <- 'precision'
  }
  else {
    structure <- covariance
    str.name <- 'covariance'
  }
  if(!is.matrix(incidence))
    stop(paste('Argument incidence must be of type matrix'))
  if(!is.matrix(structure))
    stop(paste(str.name, 'must be of type matrix'))
  if(ncol(incidence) != nrow(structure))
    stop(paste('Non conformant incidence and', str.name, 'matrices'))
  if (length(var.ini) !=1 || !is.numeric(var.ini) || var.ini <= 0)
    stop (paste('var.ini must be a positive number'))

  return(TRUE)
}

