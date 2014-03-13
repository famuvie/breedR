#' Inference with REMLF90
#' 
#' Fits a Linear Mixed Model by Restricted Maximum Likelihood
#' 
#' If either \code{genetic} and/or \code{param} are not \code{NULL}, the model 
#' residuals are assumed to have an additive genetic effects and a spatially 
#' structured random effect, respectively. The relevant parameters are 
#' \code{model} and \code{var.ini} in both cases, and \code{pedigree} in the 
#' case of a genetic effect.
#' 
#' The available models for the genetic effect are \code{add_animal}. 
#' \code{add_animal} stands for an additive animal model with a given pedigree.
#' 
#' The available models for the spatial effect are \code{Cappa07}. 
#' \code{Cappa07} uses a  two-dimensional tensor product of B-splines to 
#' represent the smooth spatially structured effect.
#' @param fixed an object of class \link{formula} (or one that can be coerced to
#'   that class): a symbolic description of the fixed effects of the model to be
#'   fitted. The details of model specification are given under 'Details'.
#' @param random if not \code{NULL}, an object of class \link{formula} with the 
#'   unstructured random effects.
#' @param genetic if not \code{NULL}, a list with relevant parameters for an 
#'   additive genetic effect; see 'Details'.
#' @param spatial if not \code{NULL}, a list with relevant parameters for a 
#'   spatial random effect; see 'Details'.
#' @param data a data frame with variables and observations
#' @param method either 'ai' or 'em' for Average-Information or 
#'   Expectation-Maximization REML respectively
#' @return An object of class 'remlf90' that can be further questioned by
#'   \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{fitted}}, etc.
#' @seealso \code{\link[pedigreemm]{pedigree}}
#' @references progsf90 wiki page: \url{http://nce.ads.uga.edu/wiki/doku.php}
#'   
#'   E. P. Cappa and R. J. C. Cantet (2007). Bayesian estimation of a surface to
#'   account for a spatial trend using penalized splines in an individual-tree 
#'   mixed model. \emph{Canadian Journal of Forest Research} 
#'   \strong{37}(12):2677-2688.
#' @export
remlf90 <- function(fixed, 
                    random  = NULL,
                    genetic = NULL, 
                    spatial = NULL,
                    data, 
                    method = c('ai', 'em')) {
  
  ## Assumptions:
  ## Only 1 pedigree
  ## Solution file header: trait effect level solution
  ## No intercept
  ## (not generalized) Linear Mixed Model
  ## There always is a genetic variance component

  ## TODO: 
  # Allow for summarized data (parameter weights)
  # Allow for multiple responses
  # Allow for generalized mixed models

  ## Checks
  if (missing(fixed) | missing(data)) { 
    stop("Usage: remlf90(fixed, data, ...); see ?remlf90")
  }
  if (class(fixed) != "formula") { 
	  stop("'fixed' should be a formula") 
  }
  if ( attr(terms(fixed), 'intercept') != 1L ) {
	  stop("There is no response in the 'fixed' argument")
  }
  if (!is.null(random)) {
	if (class(random) != "formula" | attr(terms(random), 'response') != 0L) {
		stop("random should be a response-less formula")
	}
  }
  

  #  Call
  mc <- mcout <- match.call()
  
  # Parse arguments
  method <- tolower(method)
  method <- match.arg(method)

  # Builds model frame by joining the fixed and random terms
  # and remove the intercept.
  # Add an additional 'term.types' attribute within 'terms'
  # indicating whether the term is 'fixed' or 'random'
	# progsf90 don't allow for custom model parameterizations
	# and they don't use intercepts
  mf <- build.mf(mc, remove.intercept = TRUE)
  mt <- attr(mf, 'terms')

  
  # Genetic effect
  if(!is.null(genetic)) {
    genetic$model <- match.arg(genetic$model, choices = c('add_animal'))
#     if( !inherits(genetic$pedigree, 'pedigree') )
#       stop("The argument genetic should contain a 'pedigree' object")
    if( !all(check_pedigree(genetic$pedigree)) )
      genetic$pedigree <- build_pedigree(1:3, data = genetic$pedigree)
    
    if(length(genetic$id)==1) {
      genetic$id <- data[, genetic$id]
    }
  }
  
  # Spatial effect
  if(!is.null(spatial)) {
    spatial$model <- match.arg(spatial$model,
                               choices = c('Cappa07', 'AR'))
    
    # If AR model without rho specified
    # we need to fit it with several fixed rho's
    # and return the most likely
    # TODO: It would be nice if we didn't need to recompute Q each time
    if(spatial$model == 'AR') {
      if( is.null(spatial$rho) ) spatial$rho <- matrix(c(NA, NA), 1, 2)
      if( any(is.na(spatial$rho)) ) {
        # Evaluation values for rho
        rho.grid <- build.AR.rho.grid(spatial$rho)
      } else {
        rho.grid <- spatial$rho
      }
      
      if( !is.null(nrow(spatial$rho)) ) {
        
        # Results conditional on rho
        eval.rho <- function(rho, mc) {
          mc$spatial$rho <- rho
          eval(mc)
        }
        #         test <- eval.rho(mc, c(.5, .5))
        ans.rho <- apply(rho.grid, 1, eval.rho, mc)
        # Interpolate results
        loglik.rho <- transform(rho.grid,
                                loglik = sapply(ans.rho, logLik))
        rho.idx <- which.max(loglik.rho$loglik)
        ans <- ans.rho[[rho.idx]]
        
        # Include estimation information
        ans$rho <- loglik.rho
        
        return(ans)
      }
    }
  }
  
  # Temporary files
  tmpdir <- tempdir()
  parameter.file.path <- file.path(tmpdir, 'parameters')
  if(!is.null(genetic)) genetic$tempfile <- file.path(tmpdir, 'pedigree')
  if(!is.null(spatial)) spatial$tempfile <- file.path(tmpdir, 'spatial')
  
 
  # Build a list of parameters and information for each effect
  effects <- build.effects(mf, genetic, spatial)
  
  # Generate progsf90 parameters
  # TODO: Initial variance for residuals  # FIXED ??
  # TODO: Memory efficiency. At this point there are three copies of the 
  # dataset. One in data, one in mf (only needed variables)
  # and yet one more in pf90. This is a potential problem with large datasets.
  pf90 <- progsf90(mf, effects, opt = c("sol se"), res.var.ini = 10)
  
  # Write progsf90 files
  write.progsf90(pf90, dir = tmpdir)

  # variance components and BLUPs with REML
  platform <- switch(.Platform$OS.type, 
                     unix = 'linux',
                     windows = 'windows',
                     mac = 'mac')
  binary.path <- system.file('bin', platform, package='breedR')
  
  # Change to temporal directory to avoid specification of long paths
  # Avoids Issue #1
  cdir <- setwd(tmpdir)
  reml.out <- switch(method,
                     ai = system2(file.path(binary.path, 'airemlf90'), 
                                 input = parameter.file.path, stdout=TRUE),
                     em = system2(file.path(binary.path, 'remlf90'),
                                 input = parameter.file.path, stdout = TRUE)
  )
  # Return to current directory
  setwd(cdir)
  
  
  # Error catching
  stopifnot(is.null(attr(reml.out, 'status')))

  # Parse solutions
  ans <- parse_results(file.path(tmpdir, 'solutions'), effects, mf, reml.out, method, mcout)
  
  class(ans) <- c('remlf90')  # Update to merMod in newest version of lme4 (?)
  ans
}

#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Internal methods  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Interface methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#


#' Extract fixed and random effects coefficients
#' 
#' @method coef remlf90
#' @export
coef.remlf90 <- function(object, ...) { 
  c(fixef(object), ranef(object))
}

#' Extract the Akaike Information Criterion from a fitted model
#' 
#' @S3method extractAIC remlf90
#' @export
extractAIC.remlf90 <- function(fit, scale, k,...) {
  return(fit$fit$AIC)
}
  
#' @method fitted remlf90
#' @export
fitted.remlf90 <- function (object, ...) {
  object$mu
}
  
fixef.remlf90 <- function (object, ...) {
      object$fixed
}

logLik.remlf90 <- function (object, ...) {
  # TODO: Revise this, N parameters, df, N obs.
  # I set up things such that the AIC gives what REMLF90 says
  # But I am not sure if it is the right way.
  reml.out <- object$reml$output
  rank.idx <- grep('RANK', reml.out)
  npar.idx <- grep('parameters=', reml.out)
  rank <- ifelse(identical(length(rank.idx), 1L),
                 as.numeric(strsplit(reml.out[rank.idx],
                                     split=' +RANK += +')[[1]][2]),
                 'unknown')
  npar <- ifelse(identical(length(npar.idx), 1L),
                 as.numeric(strsplit(reml.out[npar.idx],
                                     split=' # parameters= +')[[1]][2]),
                 'unknown')
  if(any(identical(rank, 'unknown') | identical(npar, 'unknown')))
    warning(paste('Could not deduce the', 
                  paste(c('rank', 'number of parameters')
                        [which(c(rank, npar)=='unknown')],
                        collapse = ' and '),
                  'from REMLF90 output'))
  
  res <- object$residual
  N <- length(res)
#   rank <- object$fit$rank
  
  if (is.null(w <- object$weights)) {
    w <- rep.int(1, N)
  }
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  ans = -object$fit[['-2logL']]/2
  attr(ans, 'df') <- npar
  attr(ans, 'nobs') <- N
  class(ans) ='logLik'
  ans
}

model.frame.remlf90 <- function (object, ...) {
  object$mf
}

model.matrtix.remlf90 <- function (object, ...) {
  object$mm
}

nobs.remlf90 <- function (object, ...) {
  nrow(as.matrix(object$y))
}

# Plotting the spatial effect in the observed locations
# TODO: implement this as a plotting method for remlf90 objects

#' Plot a model fit
#' 
#' Plots the predicted values of the spatial component of the observations' 
#' phenotypes.
#' 
#' @param x A remlf90 object, or a matrix (-like) of coordinates
#' @param y Optional. A numeric vector to be plotted.
#' @param type Character. Plot type. 'spatial' is currently the only option.
#' 
#'     
plot.remlf90 <- function (x, y = NULL, type = 'spatial', ...) {
  require(ggplot2)
  
  if( x$effects$spatial ) {
    spdat <- x$spatial$fit
    spdat$model <- x$call$spatial$model
    p <- ggplot(spdat, aes(x, y)) +
      coord_fixed() +
      geom_tile(aes(fill = z)) +
      scale_fill_gradient(low='green', high='red') +
      facet_wrap(~ model)
    p
    #       layer <- paste('geom_tile(aes(fill =', x$call$spatial$model, '))')
    #       eval(parse(text = paste('p +', layer)))
  } else {
  }
}

predict.remlf90 <- function (object, ...) {
  
}

print.remlf90 <- function (object, ...) {
  
}

#' @S3method ranef remlf90
#' @export
ranef.remlf90 <- function (object, ...) {
  object$ranef
}

residuals.remlf90 <- function (object, ...) {
  # TODO: to be used when na.action is included in remlf90
  #   naresid(object$na.action, res)
  # TODO: add a scale parameter to allow studentization
  #     first need to determine the right sigma for each residual
  object$residuals
}

#' @export
summary.remlf90 <- function(object, ...) {
  ans <- object
  
  # Literal description of the model
  effects <- paste(names(ans$effects), sep=' and ')
  title <- paste('Linear Mixed Model with', 
                 paste(effects, collapse = ' and '), 
                 ifelse(length(effects)==1, 'effect', 'effects'), 
                 'fit by', object$reml$version)
  
  # Formula
  fml <- deparse(attr(object$mf, 'terms'))
  
  # Coefficients
  # TODO: How to compute Standard errors (and therefore t scores and p-values)
  # TODO: How to avoid showing the unused levels
  coef <- do.call(rbind, object$fixed)
#   colnames(coef) <- c('Estimate')
  
  # Model fit measures
  # AIC and BIC might fail if logLik fails to retrieve
  # appropriate df or nobs attributes
  llik <- logLik(object)
  AICframe <- data.frame(AIC = tryCatch(AIC(llik), 
                                        error = function(e) object$fit$AIC), 
                         BIC = tryCatch(BIC(llik),
                                        error = function(e) 'unknown'),
                         logLik = as.vector(llik),
#TODO                          deviance = dev[["ML"]],
#                          REMLdev = dev[["REML"]],
                         row.names = "")
  ans <- c(ans, 
           model.description = title, 
           formula = fml,
           model.fit = list(AICframe),
           coefficients = list(coef)
           )
  class(ans) <- 'summary.remlf90'
  ans
}

## This is modeled a bit after  print.summary.lm :
#' @export
print.summary.remlf90 <- function(x, digits = max(3, getOption("digits") - 3),
                                  correlation = TRUE, symbolic.cor = FALSE,
                                  signif.stars = getOption("show.signif.stars"), ...) {
  
  cat(x$model.description, '\n')
  if(!is.null(x$call$formula))
    cat("Formula:", x$formula,"\n")
  if(!is.null(x$call$data))
    cat("   Data:", deparse(x$call$data),"\n")
  if(!is.null(x$call$subset))
    cat(" Subset:", x$call$subset,"\n")
  print(x$model.fit, digits = digits)
  
  if( x$effects$spatial ) {
    switch(x$spatial$name,
           AR = cat(paste("\nAutoregressive parameters for rows and columns: (",
                    paste(x$spatial$model$param, collapse = ', '),
                    ")\n", sep = '')),
           Cappa07 = cat(paste("\nNumber of inner knots for rows and columns: (",
                         paste(x$spatial$model$param, collapse =', '),
                         ")\n", sep = ''))
    )
  }
  
  cat("\nVariance components:\n")
  print(x$var, quote = FALSE, digits = digits, ...)
  if(x$effects$spatial)
    cat(" * spatial: is variance of BLUPs\n")
  
  cat('\nFixed effects:\n')
  printCoefmat(x$coefficients)
  invisible(x)
}
