#' progsf90 class
#' 
#' This function parses a model frame and extracts the relevant fields
#' that are to be written in the parameter, data and auxiliary files
#' of the progsf90 programs.
progsf90 <- function (mf, effects, opt = c("sol se"), res.var.ini = 10) {
  
  # This function parses correctly all fixed effects and the spatial effect 
  # Possibly will have to change to include competition or other complex effects
  parse.effect <- function(x) {
    if(length(x$pos) > 1){
      n.pos <- length(x$pos)
      paste(x$pos, c(rep(0, n.pos-1), n.pos), 'cov', tail(x$pos, 1) + 1:n.pos)
    }
    else  paste(do.call(c, x[c('pos', 'levels', 'type')]), collapse=' ')
  }
  
  # Build models for random effects
  random.effects.idx <- which(names(effects)=='genetic' | names(effects)=='spatial') # For the moment, the only random effects are either genetic or spatial --- TODO
  
  # Number of traits
  # (size of the response vector or matrix)
  ntraits <- ncol(as.matrix(model.response(mf)))
  
  # Weights
  weights <- ''     # No weights for the moment --- TODO
  
  # Parameters  
  par <- list(datafile = 'data',
              ntraits  = ntraits,
              neffects = sum(sapply(effects, function(x) length(x$pos))),
              observations = 1:ntraits,
              weights  = weights,
              effects  = sapply(effects, function(x) unlist(parse.effect(x))),
              residvar = res.var.ini,
              rangroup = lapply(random.effects.idx, 
                                function(x) list(pos = tail(effects[[x]]$pos, 1) - ntraits,
                                                 type = effects[[x]]$model, 
                                                 file = effects[[x]]$file, 
                                                 cov  = effects[[x]]$var)),
              options = opt
  )
  
  # Data 
  # Columns ordered as in the effects list
  # after the trait(s) 
  build.dat.single <- function(name, mf) {
    n <- nrow(mf)
    switch(name,
           '(Intercept)' = rep(1L, n),
           genetic       = effects$genetic$idx,
           spatial       = cbind(as.matrix(effects$spatial$splines$B),
                                 sapply(1:ncol(effects$spatial$splines$B),
                                        function(x) rep(x, n))),
           mf[[name]])
  }
  
  #   dat <- cbind(sapply(names(effects), build.dat.single, mf),
  #                phenotype = mf[, attr(mt, 'response')])
  dat <- do.call(cbind, 
                 c(list(phenotype = mf[, attr(attr(mf, 'terms'), 'response')]),
                   sapply(names(effects), 
                          build.dat.single, mf, simplify = FALSE)))
  
  # Additional Files
  build.file.single <- function(ef) {
    switch(ef$file,
           pedigree = list(fname = ef$file,
                          file   = ef$ped),
           spatial  = list(fname = ef$file,
                           file  = ef$splines$U))
  }
  files <- lapply(effects[random.effects.idx], build.file.single)
  
  ans <- list(parameter = par, data = dat, files = files)
  class(ans) <- 'progsf90'
  return(ans)
}


#' Build effects parameters
#' 
#' This function builds a list of effects parameters
#' as required by Misztal's progsf90 suite of programs
#' @references
#' \url{http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90.pdf}
build.effects <- function (mf, genetic, spatial) {
  
  # Build up effects data (position, levels, type)
  
  # Model terms
  mt <- attr(mf, 'terms')
  
  # Number of traits
  # (size of the response vector or matrix)
  ntraits <- ncol(as.matrix(model.response(mf)))
  
  # Position counter
  pos = ntraits + 1
  
  effects <- list()
  
  # Intercept term, if not precluded explicitly in the formula
  if(attr(mt, 'intercept')) {
    effects <- list('(Intercept)'=list(pos = pos, levels=1L, type='cross'))
    pos <- pos + 1
  }
  
  # Parameters of a single effect in the formula
  # Increases pos by one each time
  eff.par.f <- function(name) {
    # position (in the data file)
    pos <- parent.env(environment())$pos
    assign('pos', pos + 1, envir = parent.env(environment()))
    # number of levels
    nl <- ifelse(inherits(mf[[name]], 'factor'), nlevels(mf[[name]]), 1)
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
  effects <- c(effects, lapply(attr(mt, 'term.labels'), eff.par.f))
  ef_names <- attr(mt, 'term.labels')
  names(effects)[pos-ntraits - length(ef_names):1] <- ef_names
  
  
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
                        file = 'pedigree',
                        ped = as.data.frame(genetic$pedigree),
                        idx = genetic$id,
                        var = genetic$var.ini)))
    pos = pos + 1
  }
  
  # Spatial effect if applicable
  # We only have spatial coordinates of the dataset elements
  # TODO: let the user determine the degree/order of the B-splines
  
  if(!is.null(spatial)) {
    splines <- build.splines.model(spatial$coord, spatial$n.knots, degree = 3)
    effects <- c(effects, 
                 spatial = list(
                   list(pos = pos - 1 + 1:ncol(splines$B),
                        levels = 1,
                        type = 'cov',
                        model = 'user_file_i',
                        file = 'spatial',
                        var = spatial$var.ini,
                        splines = splines)))
  }
  return(effects)
}


write.progsf90 <- function (pf90, dir) {
  
  # Parameter file
  parameter.file <- 
    with(pf90$parameter, 
         c('DATAFILE', datafile, 
           'NUMBER_OF_TRAITS', ntraits,
           'NUMBER_OF_EFFECTS', neffects,
           'OBSERVATION(S)', observations,
           'WEIGHT(S)', weights,
           'EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT [EFFECT NESTED]',
           paste(unlist(effects)),
           'RANDOM_RESIDUAL VALUES', residvar,
           sapply(rangroup, 
                  function(x) c('RANDOM_GROUP', x$pos,
                                'RANDOM_TYPE', x$type,
                                'FILE', x$file,
                                '(CO)VARIANCES', x$cov)),
           paste('OPTION', options)))
  
  
  writeLines(as.character(parameter.file),
             con = file.path(dir, 'parameters'))
  # file.show(parameter.file.path)
  
  # Data file
  write.table(pf90$data,
              file = file.path(dir, pf90$parameter$datafile),
              row.names = FALSE, col.names = FALSE)
  # file.show(data.file.path)

  for(fl in pf90$files) {
    # NAs are written as 0
    write.table(fl$file, file=file.path(dir, fl$fname),
                row.names = FALSE, col.names = FALSE, na = "0")   
    # file.show(file.path(dir, fl$fname))
  }
  
  
}



#' Parse results from a progsf90 'solutions' file
parse_results <- function (solfile, effects, mf, reml.out, method, mcout) {
  
  # Parsing the results
  sol.file <- read.table(solfile, header=FALSE, skip=1)
  colnames(sol.file) <- c('trait', 'effect', 'level', 'value', 's.e.')
  
  # Assuming one trait only
  result <- by(sol.file[,4:5], sol.file$effect, identity)
  names(result) <- names(effects)
  
  # Identify factors in model terms
  mt <- attr(mf, 'terms')
  isF <- attr(mt, 'dataClasses') == 'factor' |
    attr(mt, 'dataClasses') == 'ordered'
  
  # Flags for specific effects
  isGenetic <- exists('genetic', as.environment(effects))
  isSpatial <- exists('spatial', as.environment(effects))
  
  # write labels for factor levels in results
  for( x in names(isF)[which(isF)])
    rownames(result[[x]]) <- levels(mf[[x]])
  
  # Random and Fixed effects indices
  random.effects.idx <- which(names(effects)=='genetic' | names(effects)=='spatial') # For the moment, the only random effects are either genetic or spatial --- TODO
  if( length(random.effects.idx) )
    fixed.effects.idx <- (1:length(effects))[-random.effects.idx]  else
      fixed.effects.idx <- (1:length(effects))
  
  # Fixed effects coefficients
  beta <- sol.file$value[sol.file$effect %in% fixed.effects.idx]
  #   .getXlevels(mt, mf)
  
  # Random effects coefficients
  # TODO: Return Standard Errors as well.
  # How to compute standard errors of splines predicted values?
  ranef <- list()
  if(isGenetic)
    ranef$genetic <- result$genetic$value[effects$genetic$idx]
  if(isSpatial)
    ranef$spatial <- as.vector(effects$spatial$splines$B
                               %*% 
                                 result$spatial$value)
  
  # Spatial Surface
  if (isSpatial) {
    spatial.pred <- cbind(effects$spatial$splines$plotting$grid,
                          z = as.vector(effects$spatial$splines$plotting$B
                                        %*% result$spatial$value))
  } else
    spatial.pred <- NULL
  
  # Build up the model matrix with one dummy variable per level
  # as progsf90 takes care of everything
  # I need to provide each factor with an identity matrix
  # as its 'contrasts' attribute
  diagonal_contrasts <- function(x) {
    ctr <- diag(nlevels(x))
    colnames(ctr) <- levels(x)
    attr(x, 'contrasts') <- ctr
    x
  }
  mf[isF] <- lapply(mf[isF], diagonal_contrasts)
  mm <- model.matrix(mt, mf) 
  
  eta <- drop(mm %*% beta)
  if(length(ranef)) eta <- eta + rowSums(do.call(cbind, ranef))
  
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
  
  # If spatial, report the observed spatial variance, rather than the
  # estimated variance of the spline effects which is meaningless
  if(isSpatial) { 
    varcomp['spatial', 1] <- var(ranef$spatial)
    if(method == 'ai') varcomp['spatial', 2] <- NA
  }
  
  
  # REML info
  reml.ver <- sub('^\\s+([[:graph:]]* +ver\\. +[0-9.]*).*$', '\\1', 
                  grep('REML', reml.out, value = TRUE))
  last.round.idx <- tail(grep('In round', reml.out), 1)
  last.round <- as.numeric(strsplit(strsplit(reml.out[last.round.idx],
                                             split='In round')[[1]][2],
                                    split='convergence=')[[1]])
  reml <- list(
    version = gsub('\\s+', ' ', reml.ver),
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
  y <- model.response(mf, "numeric")
  
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
    effects = list(pedigree = isGenetic,
                   spatial  = isSpatial), # TODO competition, ...
    mf = mf,
    mm = mm,
    y = y,
    fixed = result[fixed.effects.idx],
    ranef = ranef,
    eta = eta,
    mu = mu,
    residuals = y - mu,
    spatial = list(model      = effects$spatial$splines[1:3],
                   prediction = spatial.pred),
    var = varcomp,
    fit = fit,
    reml = reml
  )
  
}
