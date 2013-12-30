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
#'    \code{\link{http://nce.ads.uga.edu/wiki/doku.php}}
#'    
#'    E. P. Cappa and R. J. C. Cantet (2007). Bayesian estimation of a surface to account for a spatial trend using penalized splines in an individual-tree mixed model. \emph{Canadian Journal of Forest Research} \strong{37}(12):2677-2688.
#' @export
remlf90 <- function(formula, genetic=NULL, spatial=NULL, data, id = 1:nrow(data), method=c('ai', 'em')) {
  
  ## Assumptions:
  ## Only 1 pedigree
  ## Solution file header: trait effect level solution
  ## No intercept
  ## (not generalized) Linear Mixed Model
  ## There always is a genetic variance component
  
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
  
  if(length(id)==1) id <- data[, id]

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
    ped <- as.data.frame(genetic$pedigree)   # Extract pedigree
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
  
  # Build models for random effects
  random.effects.idx <- which(names(effects)=='genetic' | names(effects)=='spatial') # For the moment, the only random effects are either genetic or spatial --- TODO
  
  
  
  # Number of traits
  # (size of the response vector or matrix)
  ntraits <- ncol(as.matrix(model.response(mf)))
  
  # Write the parameter file
  weights <- ''     # No weights for the moment --- TODO
  res.var.ini <- 10 # Initial variance for residuals  # FIXED ??
  parameter.file <- c('DATAFILE', data.file.path, 
                      'NUMBER_OF_TRAITS', ntraits, 
                      'NUMBER_OF_EFFECTS', length(effects),
                      'OBSERVATION(S)', length(effects) + 1:ntraits,
                      'WEIGHT(S)', weights,
                      'EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT [EFFECT NESTED]',
                      paste(lapply(effects, function(x) paste(do.call(c, x[c('pos', 'levels', 'type')]), collapse=' '))),
                      'RANDOM_RESIDUAL VALUES', res.var.ini, 
                      sapply(random.effects.idx, function(x) c('RANDOM_GROUP', x, 'RANDOM_TYPE', effects[[x]]$model, 'FILE', effects[[x]]$file, '(CO)VARIANCES', effects[[x]]$var))
  )
  
  writeLines(parameter.file, con = parameter.file.path)
  # file.show(parameter.file.path)
  
  # Write the data file
  # Columns ordered as in the effects list
  # TODO data[, id] in build.dat.single will fail if the data is not provided
  # everything we need we have to take it from mf.
  build.dat.single <- function(name, mf) {
    n <- nrow(mf)
    switch(name,
           '(Intercept)' = rep(1L, n),
           genetic       = id,
           spatial       = 1:n,
           mf[[name]])
  }
  
  dat <- cbind(sapply(names(effects), build.dat.single, mf),
               phenotype = mf[, attr(mt, 'response')])
  
  write.table(dat, file = data.file.path, row.names = FALSE, col.names = FALSE)
  # file.show(data.file.path)
  
  # Write the pedigree file if applicable
  if(!is.null(genetic)) {
    # ASSUMPTION: there is only one effect with an associated pedigree
    write.table(ped, file=genetic$tempfile, row.names = FALSE, col.names = FALSE, na = "0")   # NAs are written as 0
    # file.show(pedigree.file.path)
  }
  
  # Write the spatial file if applicable
  if(!is.null(spatial)) {
    # Write U matrix
  }
  
  # variance components with REML
  platform <- switch(.Platform$OS.type, 
                     unix = 'linux',
                     windows = 'windows',
                     mac = 'mac')
  binary.path <- system.file('bin', platform, package='breedR')
  reml.out <- switch(method,
                     ai = system(shQuote(file.path(binary.path, 'airemlf90')), input = parameter.file.path, intern=TRUE),
                     em = system(shQuote(file.path(binary.path, 'remlf90')), input = parameter.file.path, intern=TRUE)
  )
  
  # Error catching
  stopifnot(is.null(attr(reml.out, 'status')))
  
  # Parsing the results
  sol.file <- read.table('solutions', header=FALSE, skip=1)
  colnames(sol.file) <- c('trait', 'effect', 'level', 'value')
  file.remove('solutions')
  
  # One trait only
  result <- tapply(sol.file$value, sol.file$effect, identity)
  names(result) <- names(effects)

  
  # Fixed effects coefficients
  fixed.effects.idx <- (1:length(effects))[-random.effects.idx]
  beta <- sol.file$value[sol.file$effect %in% fixed.effects.idx]
#   .getXlevels(mt, mf)
  
  
  # INTERCEPT ISSUE
  # PROGSF90 do not use an intercept when it has factor variables
  # instead, they fit a parameter for every level of every factor.
  # Therefore, we compute the model matrix with the intercept removed.
  mm <- model.matrix(update(formula, ~ . -1), mf)
  # Alternatively, we introduce an intercept term in the reml data
  # But as a result, the random effects remain the same, and the
  # fixed effects take the *last* level as the reference (effect=0)
  mm <- model.matrix(formula, mf)
  
  # Linear Predictor
  # ASSUMPTION: the model has no intercept
  # TODO: This can be computed more efficiently by using the model matrix
  # as soon as I can match the progsf90 parameterization with the 
  # standard parameterization in R
  eta.i <- function(levels) {
    sum(mapply(function(y, x) y[x], result, levels))
  }
  eta = apply(dat[,-(length(effects) + 1:ntraits)], 1, eta.i)
  
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
  
  gen.var.idx <- grep('Genetic variance', reml.out) + 1
  stopifnot(length(gen.var.idx)==1)
  gen.var <- as.numeric(reml.out[gen.var.idx])

  gen.var.sd.idx <- grep(paste(sd.label, 'for G'), reml.out) + 1
  stopifnot(length(gen.var.sd.idx)==1)
  gen.var.sd <- as.numeric(reml.out[gen.var.sd.idx])
  
  res.var.idx <- grep('Residual variance', reml.out) + 1
  stopifnot(length(res.var.idx)==1)
  res.var <- as.numeric(reml.out[res.var.idx])
  
  res.var.sd.idx <- grep(paste(sd.label, 'for R'), reml.out) + 1
  stopifnot(length(res.var.sd.idx)==1)
  res.var.sd <- as.numeric(reml.out[res.var.sd.idx])
  
  
  # REML info
  last.round.idx <- tail(grep('In round', reml.out), 1)
  last.round <- as.numeric(strsplit(strsplit(reml.out[last.round.idx], split='In round')[[1]][2], split='convergence=')[[1]])
  reml <- list(
    rounds = last.round[1],
    convergence = last.round[2],
    delta.conv = as.numeric(strsplit(reml.out[last.round.idx+1], split='delta convergence=')[[1]][2]),
    output = reml.out
    )
  
  # Fit info
  last.fit <- as.numeric(strsplit(strsplit(reml.out[last.round.idx-1], split='-2logL =')[[1]][2], split=': AIC =')[[1]])
  fit <- list(
    '-2logL' = last.fit[1],
    AIC = last.fit[2]
    )
  
  # Observed response
  y = dat[, length(effects) + 1:ntraits]
  
#   # Response in the linear predictor scale
#   # TODO: apply link
#   y.scaled <- y
  
  ans <- list(
    call = mc,
    method = method,
    effects = list(pedigree = !is.null(genetic)), # TODO spatial, competition, ...
    mf = mf,
    mm = 'TODO',
    y = y,
    fixed = result[-random.effects.idx],
    ranef = result[random.effects.idx],
    eta = eta,
    mu = mu,
    residuals = y - mu,
    var = rbind(genetic=c(mean=gen.var, sd=gen.var.sd), 
                residual=c(res.var, res.var.sd)),
    fit = fit,
    reml = reml
    )
  class(ans) <- c('remlf90')  # Update to merMod in newest version of lme4 (?)
  ans
}

#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Internal methods  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#

#' Build effects parameters
#' 
#' This function builds a list of effects parameters
#' as required by Mistal'z progsf90 suite of programs
#' @references
#' \url{http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90.pdf}
build.effects <- function (mf, genetic, spatial) {
  
  # Build up effects data (position, levels, type)
  
  # Model terms
  mt <- attr(mf, 'terms')

  # Intercept term, if not precluded explicitly in the formula
  if(attr(mt, 'intercept')) effects <- list('(Intercept)'=list(pos = 1L, levels=1L, type='cross'))
  
  # Parameters of a single effect in the formula
  eff.par.f <- function(name) {
    # position (in the data file)
    pos <- parent.env(environment())$pos
    assign('pos', pos + 1, envir = parent.env(environment()))
    # number of levels
    nl <- ifelse(inherits(mf[[name]], 'factor'), nlevels(mf[[name]]), length(unique(mf[[name]])))
    # type: factors = "cross"; continuous = "cov"
    type <- switch(attr(mt, 'dataClasses')[name],
                   ordered = 'cross',
                   factor = 'cross',
                   numeric = 'cov',
                   'cross')
    # nested effects
    # TODO
    return(list(pos=pos, levels=nl, type=type))
  }
  
  # Parameters for all effects in the formula
  pos = 2L
  effects <- c(effects, lapply(attr(mt, 'term.labels'), eff.par.f))
  names(effects)[-1] <- attr(mt, 'term.labels')
  
  # Genetic effect
  # Both the genetic and spatial terms are "cross"
  # For pedigree effects, there might be more levels than those
  # present in the data. We should declare the levels present in the pedigree.
  if(!is.null(genetic)) {
    effects <- c(effects, 
                 genetic = list(
                   list(pos = pos,
                        levels = nrow(as.data.frame(genetic$pedigree)),
                        type = 'cross',
                        model = genetic$model,
                        file = genetic$tempfile,
                        var = genetic$var.ini)))
    pos = pos + 1
  }
  
  # Spatial effect
  # We only have spatial coordinates of the dataset elements
  if(!is.null(spatial)) {
    effects <- c(effects, 
                 spatial = list(
                   list(pos = pos,
                        levels = nrow(mf),
                        type = 'cross',
                        model = 'user_file',
                        file = spatial$tempfile,
                        var = spatial$var.ini)))
  }
  return(effects)
}



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
extractAIC.remlf90 <- function(object, scale, k,...) {
  
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
  rank <- as.numeric(strsplit(reml.out[rank.idx], split=' +RANK += +')[[1]][2])
  npar.idx <- grep('parameters=', reml.out)
  npar <- as.numeric(strsplit(reml.out[npar.idx], split=' # parameters= +')[[1]][2])
  
  res <- object$residual
  N <- length(res)
  rank <- object$fit$rank
  
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
  title <- paste('Linear Mixed Model with', effects, ifelse(length(effects)==1, 'effect', 'effects'), 'fit by', paste(toupper(ans$method), 'REML', sep='-'))
  
  # Formula
  fml <- deparse(attr(object$mf, 'terms'))
  
  # Coefficients
  # TODO: How to compute Standard errors (and therefore t scores and p-values)
  # TODO: How to avoid showing the unused levels
  coef <- as.matrix(unlist(object$fixed))
  colnames(coef) <- c('Estimate')
  
  # Model fit measures
  llik <- logLik(object)
  AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
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
