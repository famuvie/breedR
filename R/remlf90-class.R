#### Inference with REMLF90 ####
remlf90 <- function(formula, effects, data, method=c('ai', 'em')) {
  
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
  
  # TODO: Allow for subsetting the data (i.e. exclude founders)
  # in the function call (parameter subset). Use this parameter
  # when building up the model frame.
  
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
  # and the formula must be evaluated in the calling environment (parent)
  # besides, it allows passing additional arguments to model.frame
  # like subset, na.action, etc.
  mc[[1]] <- quote(stats::model.frame)
  mc$effects <- NULL
  mf <- eval(mc, parent.frame(1L))
  mt <- attr(mf, 'terms')
#   mf <- model.frame(update(formula, ~.-1), data)
  # Better add an intercept to progsf90
  
#   # Extract data columns in the right order
#   # (the same as in the formula)
#   # the response goes to the last column
#   tf <- terms.formula(formula)
#   dat <- data[,c(attr(attr(mf, 'terms'), 'term.labels'), as.character(tf[[2]]))]
#   # This is bullshit. Do it right with model.frame and so on
#   # TODO: adapt everything

  # Extract pedigree
  # TODO: If there is. Include a parameter pedigree
  withPedigree = TRUE
  ped.effect.idx <- which(sapply(effects, function(x) !is.null(x$pedigree)))
  ped <- as.data.frame(effects[[ped.effect.idx]]$pedigree)
  
  # Temporary files
  tmpdir <- tempdir()
  parameter.file.path <- file.path(tmpdir, 'parameters')
  data.file.path <- file.path(tmpdir, 'data')
  pedigree.file.path <- file.path(tmpdir, 'pedigree')
  
  # Complete effects data (position and levels)
  # The order in the list of effects is *assumed*
  # the same as in the formula terms

  # the given effects list must have the same number of therms as the formula
  stopifnot(identical(length(effects), length(attr(mt, 'term.labels'))))
  
  
  # Levels of factors
  # Watch out! the parameters file must read the levels present in the pedigree
  # even if the data is restricted
  factor.effects <- which(sapply(effects, function(x) x$type == 'cross'))
    # TODO: We should better assess factors from data (model frame)
    # rather than from this artificial effects argument
    # TODO: We also need the levels of numeric (cov) covariates
  factor.levels <- sapply(as.data.frame(as.matrix(mf[,-1])[,factor.effects]), function(x) nlevels(as.factor(x)))
    # NOTE: as.data.frame and as.matrix are needed to avoid simplification in the case there is only one effect.
  effects[] <- mapply(c, lapply(factor.levels, function(x) list(levels=x)), effects, SIMPLIFY=FALSE)
  # For pedigree effects, there might be more levels than those
  # present in the data. We should state the levels present in the pedigree.
  effects[[ped.effect.idx]]$levels <- nrow(ped)
  # Add intercept to effects list, if not precluded explicitly in the formula
  # Caution: this must be after adding the levels and befor adding the position
  if(attr(mt, 'intercept')) effects <- c(list('(Intercept)'=list(levels=1L, type='cross', model='fixed')), effects)
  # Position in formula
  effects[] <- mapply(c, lapply(1:length(effects), function(x) list(pos=x)), effects)
#   # Retrieve lost names # Unnecessary?
#   effects[] <- ... preserves the names!
  
  
  # Number of traits
  # (size of the response vector)
  # TODO: this is probably wrong. If the response was a matrix,
  # there would still be one single term in the formula.
  # Better use (check later)
#   ncol(as.matrix(model.response(mf)))
  ntraits <- length(formula[[2]])
  
  # Write the parameter file
  weights <- ''     # No weights for the moment
  res.var.ini <- 10 # Initial variance for residuals  # FIXED ??
  random.effects.idx <- which(sapply(effects, function(x) x$model != 'fixed'))
  parameter.file <- c('DATAFILE', data.file.path, 
                      'NUMBER_OF_TRAITS', ntraits, 
                      'NUMBER_OF_EFFECTS', length(effects),
                      'OBSERVATION(S)', length(effects) + 1:ntraits,
                      'WEIGHT(S)', weights,
                      'EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT [EFFECT NESTED]',
                      paste(lapply(effects, function(x) paste(do.call(c, x[c('pos', 'levels', 'type')]), collapse=' '))),
                      'RANDOM_RESIDUAL VALUES', res.var.ini, 
                      sapply(random.effects.idx, function(x) c('RANDOM_GROUP', x, 'RANDOM_TYPE', effects[[x]]$model, 'FILE', pedigree.file.path, '(CO)VARIANCES', effects[[x]]$var))
  )
  
  writeLines(parameter.file, con = parameter.file.path)
  # file.show(parameter.file.path)
  
  # Write the data file
  # Columns ordered as in the effects list (or formula)
  # TODO: derive this from the formula
  # strsplit(labels(terms(formula)), '^.*\\(')
  # match(labels(terms(formula)), names(data(meta)))
  # sapply(names(data(meta)), grep, labels(terms(formula)))
  # In the simplest case, the terms of the formula correspond to 
  # variable names in the data file
  # dat <- head(dat)
  # Numerically encode factors
  dat.factor.idx <- sapply(mf, is.factor)
  dat <- mf
  dat[dat.factor.idx] <- sapply(mf[dat.factor.idx], unclass)
  # Add the intercept column, and put the response(s) last
  dat <- cbind(data.frame(intercept=1L), dat[, -attr(mt, 'response')], dat[, attr(mt, 'response')])
  
  write.table(dat, file = data.file.path, row.names = FALSE, col.names = FALSE)
  # file.show(data.file.path)
  
  # Write the pedigree file
  # ASSUMPTION: there is only one effect with an associated pedigree
  write.table(ped, file=pedigree.file.path, row.names = FALSE, col.names = FALSE, na = "0")   # NAs are written as 0
  # file.show(pedigree.file.path)
  
  # variance components with REML
  reml.out <- switch(method,
                     ai = system('~/t4f/bin/airemlf90', input = parameter.file.path, intern=TRUE),
                     em = system('~/t4f/bin/remlf90', input = parameter.file.path, intern=TRUE)
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
  gen.var.idx <- grep('Genetic variance', reml.out) + 1
  stopifnot(length(gen.var.idx)==1)
  gen.var <- as.numeric(reml.out[gen.var.idx])

  gen.var.sd.idx <- grep('S.D. for G', reml.out) + 1
  stopifnot(length(gen.var.sd.idx)==1)
  gen.var.sd <- as.numeric(reml.out[gen.var.sd.idx])
  
  res.var.idx <- grep('Residual variance', reml.out) + 1
  stopifnot(length(res.var.idx)==1)
  res.var <- as.numeric(reml.out[res.var.idx])
  
  res.var.sd.idx <- grep('S.D. for R', reml.out) + 1
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
    effects = list(pedigree = withPedigree), # TODO spatial, competition, ...
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

#### Non-exported methods ####

coef.remlf90 <- function(object, ...) { 
  c(fixef(object), ranef(object))
}

extractAIC.remlf90 <- function(object, ...) {
  
}
  
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
  res.f90$ranef
}

residuals.remlf90 <- function (object, ...) {
  # TODO: to be used when na.action is included in remlf90
  #   naresid(object$na.action, res)
  # TODO: add a scale parameter to allow studentization
  #     first need to determine the right sigma for each residual
  object$residuals
}

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
print.summary.remlf90 <- function(x, digits = max(3, getOption("digits") - 3),
                                  correlation = TRUE, symbolic.cor = FALSE,
                                  signif.stars = getOption("show.signif.stars"), ...) {
  llik <- x$fit$'-2logL'
  
  cat(x$model.description, '\n')
  if(!is.null(x$call$formula))
    cat("Formula:", x$formula,"\n")
  if(!is.null(x$call$data))
    cat("   Data:", x$call$data,"\n")
  if(!is.null(x$call$subset))
    cat(" Subset:", x$call$subset,"\n")
  print(x$model.fit, digits = digits)
  
  cat("\nVariance components:\n")
  print(x$var, quote = FALSE, digits = digits, ...)
  
  cat('\nFixed effects:\n')
  printCoefmat(x$coefficients)
  invisible(x)
}
