#' Inference with REMLF90
#' 
#' Fits a Linear Mixed Model by Restricted Maximum Likelihood
#' @param formula an object of class \link{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param genetic if not \code{NULL}, a list with relevant parameters for an additive genetic effect; see 'Details'.
#' @param spatial if not \code{NULL}, a list with relevant parameters for a spatial random effect; see 'Details'.
#' @details 
#' If either \code{genetic} and/or \code{param} are not \code{NULL}, the model residuals are assumed to have an additive genetic effects and a spatially structured random effect, respectively.
#' The relevant parameters are \code{model} and \code{var.ini} in both cases, and \code{pedigree} in the case of a genetic effect.
#' 
#' The available models for the genetic effect are \code{add_animal}. \code{add_animal} stands for an additive animal model with a given pedigree.
#' The available models for the spatial effect are \code{Cappa07}. \code{Cappa07} uses a  two-dimensional tensor product of B-splines to represent the smooth spatially structured effect.
#' @seealso \code{\link[pedigreemm]{pedigree}}
#' @references
#'    progsf90 wiki page: \url{http://nce.ads.uga.edu/wiki/doku.php}
#'    
#'    E. P. Cappa and R. J. C. Cantet (2007). Bayesian estimation of a surface to account for a spatial trend using penalized splines in an individual-tree mixed model. \emph{Canadian Journal of Forest Research} \strong{37}(12):2677-2688.
#' @export
remlf90 <- function(formula, genetic=NULL, spatial=NULL, data, method=c('ai', 'em')) {
  
  ## Assumptions:
  ## Only 1 pedigree
  ## Solution file header: trait effect level solution
  ## No intercept
  ## (not generalized) Linear Mixed Model
  ## There always is a genetic variance component

  # TODO: Allow for removing the intercept in the formula
  
  # TODO: Allow for other (diagonal) random effects (notation??)
  
  # TODO: Allow for summarized data (parameter weights)
  
  # TODO: discriminate properly factor covariates (cross) 
  # from numeric covariates (cov)
  
  # TODO: Allow for multiple responses
  
  # TODO: Allow for generalized mixed models
  
  #  Call
  mc <- mcout <- match.call()
  
  # Parse arguments
  method <- tolower(method)
  method <- match.arg(method)

  ## parse data and formula
  # NOTE: This complicated way of calling a function accounts
  # for the fact that the user possibly didn't pass a data argument
  # and the formula must be evaluated in the calling environment (parent.frame())
  # besides, it allows passing additional arguments to model.frame
  # like subset, na.action, etc.
  mc[[1]] <- quote(stats::model.frame)
  mc$genetic <- mc$spatial <- mc$method <- NULL
  mf <- eval(mc, parent.frame())
  mt <- attr(mf, 'terms')
#   mf <- model.frame(update(formula, ~.-1), data)
  # Better add an intercept to progsf90
  
  # Genetic effect
  if(!is.null(genetic)) {
    genetic$model <- match.arg(genetic$model, choices = c('add_animal'))
    if(!inherits(genetic$pedigree, 'pedigree')) stop("The argument genetic should contain a 'pedigree' object")
    if(length(genetic$id)==1) {
      # TODO: Do it right. data need not be present.
#       mc[[2]] <- ~ genetic$id
#       genetic$id <- eval(mc, parent.frame())[[1]]
      genetic$id <- data[, genetic$id]
    }
  }
  
  # Spatial effect
  if(!is.null(spatial)) {
    spatial$model <- match.arg(spatial$model, choices = c('Cappa07'))
  }
  
  # Temporary files
  tmpdir <- tempdir()
  # WORKAROUND: in Windows, tmpdir is too lengthy for AIREMLF90 
  # This is somewhat dangerous: there might be permission issues
  # Fixes Issue #1
  if(.Platform$OS.type == 'windows') {
    tmpdir = "C:\\Rtmp"
    if(file.exists(tmpdir)) stop(paste(tmpdir, 'already exists'))
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive=TRUE))
  }
  parameter.file.path <- file.path(tmpdir, 'parameters')
  data.file.path <- file.path(tmpdir, 'data')
  if(!is.null(genetic)) genetic$tempfile <- file.path(tmpdir, 'pedigree')
  if(!is.null(spatial)) spatial$tempfile <- file.path(tmpdir, 'spatial')
  
 
  # Build effects' parameters
  effects <- build.effects(mf, genetic, spatial)
  
  # Generate progsf90 parameters
  # Initial variance for residuals  # FIXED ??
  pf90 <- progsf90(mf, effects, res.var.ini = 10)
  
  # Write progsf90 files
  write.progsf90(pf90, dir = tmpdir)


  
  # variance components with REML
  platform <- switch(.Platform$OS.type, 
                     unix = 'linux',
                     windows = 'windows',
                     mac = 'mac')
  binary.path <- system.file('bin', platform, package='breedR')
  
  # Change to temporal directory to avoid specification of long paths
  # Avoids Issue #1
  cdir <- setwd(tmpdir)
  reml.out <- switch(method,
                     ai = system(shQuote(file.path(binary.path, 'airemlf90')), input = parameter.file.path, intern=TRUE),
                     em = system(shQuote(file.path(binary.path, 'remlf90')), input = parameter.file.path, intern=TRUE)
  )
  # Return to current directory
  setwd(cdir)
  
  
  # Error catching
  stopifnot(is.null(attr(reml.out, 'status')))
  
  # Parsing the results
  sol.file <- read.table(file.path(tmpdir, 'solutions'), header=FALSE, skip=1)
  colnames(sol.file) <- c('trait', 'effect', 'level', 'value', 's.e.')
  
  # One trait only
  result <- by(sol.file[,4:5], sol.file$effect, identity)
  names(result) <- names(effects)

  
  # Random and Fixed effects indices
  random.effects.idx <- which(names(effects)=='genetic' | names(effects)=='spatial') # For the moment, the only random effects are either genetic or spatial --- TODO
  fixed.effects.idx <- (1:length(effects))[-random.effects.idx]
  
  # Fixed effects coefficients
  beta <- sol.file$value[sol.file$effect %in% fixed.effects.idx]
#   .getXlevels(mt, mf)
  
  # Random effects coefficients
  # TODO: Return Standard Errors as well.
  # How to compute standard errors of splines predicted values?
  ranef <- list()
  if(!is.null(genetic))
    ranef$genetic <- result$genetic$value[genetic$id]
  if(!is.null(spatial))
    ranef$spatial <- as.vector(effects$spatial$splines$B
                               %*% 
                               result$spatial$value)
  
  # Spatial Surface
  if (!is.null(spatial)) {
    spatial.pred <- cbind(effects$spatial$splines$plotting$grid,
                          z = as.vector(effects$spatial$splines$plotting$B
                                        %*% result$spatial$value))
  } else
    spatial.pred <- NULL
  
  # INTERCEPT ISSUE
  # PROGSF90 do not use an intercept when it has factor variables
  # instead, they fit a parameter for every level of every factor.
  # Therefore, we compute the model matrix with the intercept removed.
  mm <- model.matrix(update(formula, ~ . -1), mf)
  # Alternatively, we introduce an intercept term in the reml data
  # But as a result, the random effects remain the same, and the
  # fixed effects take the *last* level as the reference (effect=0)
  # Watch out! Sometimes Misztal takes the intercept as reference
  # I don't know the rule for this.
#   mm <- model.matrix(formula, mf)
  
#   # Linear Predictor
#   # ASSUMPTION: the model has no intercept
#   # TODO: This can be computed more efficiently by using the model matrix
#   # as soon as I can match the progsf90 parameterization with the 
#   # standard parameterization in R
#   eta.i <- function(levels) {
#     sum(mapply(function(y, x) y[x], result, levels))
#   }
#   eta = apply(dat[,-(length(effects) + 1:ntraits)], 1, eta.i)

  # TODO: Misztal uses either the intercept or one of the levels
  # of a factor as reference (so one of the elements in beta is zero)
  # But I don't know which one it will be.
  # With globulus it is the intercept, while with m4 is the gen4.
  # For the moment, I build manually a full model matrix
  mm <- cbind(1, mm)
#   # Exactly one level should be zeroed
#   # Update: Not true. Omitting this check.
#   stopifnot(identical(sum(sapply(beta, identical, 0)), 1L))
  eta.genetic <- eta.spatial <- 0
  eta <- mm %*% beta + rowSums(do.call(cbind, ranef))
  
  # Fitted Values
  # ASSUMPTION: Linear Model (not generalized)
  # TODO: apply inverse link
  mu = eta
  
  # Variance components
  # ASSUMPTION: I always have a Genetic Variance component
  
  # Issue #2 In Linux, AIREMLF90 prints S.D. for R and G
  # while under Windows it outputs SE for R and G
  # Update: from version 1.109 (at least), Linux updated to SE as well
#   sd.label <- ifelse(.Platform$OS.type == 'windows', 'SE', 'S.D.')
  sd.label <- ifelse(TRUE, 'SE', 'S.D.')
  
  varcomp.idx <- grep('Genetic variance|Residual variance', reml.out) + 1
  # There should be one variance for each random effect plus one resid. var.
  stopifnot(identical(length(varcomp.idx), length(random.effects.idx) + 1L))
  varcomp <- as.numeric(reml.out[varcomp.idx])
  names(varcomp) <- c(names(effects)[random.effects.idx], 'residual')
  varcomp <- cbind('Estimated variances' = varcomp)
  
  # REML does not print Standard Errors for variance components
  if(method == 'ai'){
    varsd.idx <- grep(paste(sd.label, 'for G|for R'), reml.out) + 1
    # There should be one variance for each random effect plus one resid. var.
    stopifnot(identical(length(varcomp.idx), length(random.effects.idx) + 1L))
    varcomp <- cbind(varcomp, 'S.E.' = as.numeric(reml.out[varsd.idx]))
  }
  

  
  # REML info
  reml.ver <- grep('REML', reml.out, value = TRUE)
  last.round.idx <- tail(grep('In round', reml.out), 1)
  last.round <- as.numeric(strsplit(strsplit(reml.out[last.round.idx],
                                             split='In round')[[1]][2],
                                    split='convergence=')[[1]])
  reml <- list(
    version = reml.ver,
    rounds = last.round[1],
    convergence = last.round[2],
    delta.conv = as.numeric(strsplit(reml.out[last.round.idx+1],
                                     split='delta convergence=')[[1]][2]),
    output = reml.out
    )
  
  # Fit info
  last.fit <- as.numeric(strsplit(strsplit(reml.out[last.round.idx-1],
                                           split='-2logL =')[[1]][2],
                                  split=': AIC =')[[1]])
  fit <- list(
    '-2logL' = last.fit[1],
    AIC = last.fit[2]
    )
  
  # Observed response
  y = mf[[attr(attr(mf, 'terms'), 'response')]]
  
#   # Response in the linear predictor scale
#   # TODO: apply link
#   y.scaled <- y
  
# TODO: Add inbreeding coefficient for each tree (in the pedigree) (use pedigreemm::inbreeding())
#       Include the matrix A of additive relationships (as sparse. Use pedigreemm::getA)
#       Compute the heritability estimates and its standard error 
#       Compute covariances estimates for multiple traits (and their standard errors)
  
  ans <- list(
    call = mcout,
    method = method,
    effects = list(pedigree = !is.null(genetic),
                   spatial  = !is.null(spatial)), # TODO competition, ...
    mf = mf,
    mm = mm,
    y = y,
    fixed = result[fixed.effects.idx],
    ranef = ranef,
    eta = eta,
    mu = mu,
    residuals = y - mu,
    spatial = spatial.pred,
    var = varcomp,
    fit = fit,
    reml = reml
    )
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
#' @method extractAIC remlf90
#' @export
extractAIC.remlf90 <- function(fit, scale, k,...) {
  
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
  N0 <- N
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

plot.remlf90 <- function (object, ...) {
  
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
                 'fit by', paste(toupper(ans$method), 'REML', sep='-'))
  
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
                                        error = function(e) 'unknown'), 
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
  llik <- x$fit$'-2logL'
  
  cat(x$model.description, '\n')
  if(!is.null(x$call$formula))
    cat("Formula:", x$formula,"\n")
  if(!is.null(x$call$data))
    cat("   Data:", deparse(x$call$data),"\n")
  if(!is.null(x$call$subset))
    cat(" Subset:", x$call$subset,"\n")
  print(x$model.fit, digits = digits)
  
  cat("\nVariance components:\n")
  print(x$var, quote = FALSE, digits = digits, ...)
  
  cat('\nFixed effects:\n')
  printCoefmat(x$coefficients)
  invisible(x)
}
