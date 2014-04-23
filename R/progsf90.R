#' progsf90 class
#' 
#' This function parses a model frame and extracts the relevant fields
#' that are to be written in the parameter, data and auxiliary files
#' of the progsf90 programs.
progsf90 <- function (mf, effects, opt = c("sol se"), res.var.ini = 10) {
  
  # Parses all fixed; random; spatial and genetic effects
  # Possibly it will need changes to account for competition 
  # or other complex effects.
  # Builds the lines in the EFFECTS section
  parse.effect <- function(x) {
    if(length(x$pos) > 1) {
      # Effects that reflect into multiple lines (e.g. spatial)
      n.pos <- length(x$pos)
      paste(x$pos, c(rep(0, n.pos-1), n.pos), 'cov', tail(x$pos, 1) + 1:n.pos)
    }
    else {
      # Simple effects
      paste(do.call(c, x[c('pos', 'levels', 'type')]), collapse=' ')
    }
  }
  
  # Build models for random effects
  mt <- attr(mf, 'terms')
  random.effects.idx <- c(which(attr(mt, 'term.types') == 'random'),
                          which(names(effects) == 'genetic' | 
                                  names(effects)=='spatial'))
  
                          
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
              effects  = sapply(effects,
                                function(x) unlist(parse.effect(x))),
              residvar = res.var.ini,
              rangroup = lapply(random.effects.idx, 
                                function(x) {
                                  list(pos = tail(effects[[x]]$pos, 1)-ntraits,
                                       type = effects[[x]]$model, 
                                       file = effects[[x]]$file, 
                                       cov  = effects[[x]]$var)
                                }),
              options = opt
  )
  
  # Data 
  # Columns ordered as in the effects list
  # after the trait(s) 
  build.dat.single <- function(name, mf) {
    n <- nrow(mf)
    # Distinguish between the splines model and the AR model
    process.spatial.dat <- function(x) {
      if(length(x$pos) == 1) { # AR
        x$sp$B
      } else {
        # For the splines, I need both the columns of B
        # and further columns like
        # 1  2  3  ...
        # 1  2  3  ...
        # ·  ·  ·  ...
        # ·  ·  ·  ...
        cbind(as.matrix(x$sp$B),
              sapply(1:ncol(x$sp$B),
                     function(y) rep(y, n)))
      }
    }
    switch(name,
           '(Intercept)' = rep(1L, n),
           genetic       = effects$genetic$idx,
           spatial       = process.spatial.dat(effects$spatial),
           mf[[name]])
  }
  
  # Phenotype. Account for missing values.
  # Encode NA as 0 for Misztal's programs. Make sure there are no "real" zeroes.
  Y <- mf[, attr(attr(mf, 'terms'), 'response')]
  if( any(is.na(Y)) ) {
    if( 0. %in% Y )
      stop("Can't include missing values while some real observations are 0.\n")
    Y[which(is.na(Y))] <- 0
  }
  
  dat <- do.call(cbind, 
                 c(list(phenotype = Y),
                   sapply(names(effects), 
                          build.dat.single, mf, simplify = FALSE)))
  
  # Additional Files
  build.file.single <- function(ef) {
    switch(ef$file,
           pedigree = list(fname = ef$file,
                          file   = ef$ped),
           spatial  = list(fname = ef$file,
                           file  = ef$sp$U),
           NULL)
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
build.effects <- function (mf, genetic, spatial, var.ini) {
  
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
  
  # Parameters for every single effect in the formula
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
  
  # Unstructured random effects
  # We need to specify 'model', 'file' and 'var'
  rnd.idx <- which(attr(mt, 'term.types') == 'random')
  for (term in rnd.idx) {
    effects[[term]] <- c(effects[[term]], 
                         list(model = 'diagonal',
                              file  = '',
                              var   = var.ini[[names(effects)[term]]]))
  }
  
  # Genetic effect
  # Both the genetic and spatial terms are "cross"
  # For pedigree effects, there might be more levels than those
  # present in the data. We should declare the levels present in the pedigree.
  # It is possible that the pedigree has been recoded/reordered
  # In that case, we need to recode the data file id codes as well
  if( !is.null(attr(genetic$pedigree, 'map')) )
    idx <- attr(genetic$pedigree, 'map')[genetic$id]
  else
    idx <- genetic$id
  
  if(!is.null(genetic)) {
    effects <- c(effects, 
                 genetic = list(
                   list(pos = pos,
                        levels = nrow(as.data.frame(genetic$pedigree)),
                        type = 'cross',
                        model = genetic$model,
                        file = 'pedigree',
                        ped = as.data.frame(genetic$pedigree),
                        idx = idx,
                        var = genetic$var.ini)))
    pos = pos + 1
  }
  
  # Spatial effect if applicable
  # We only have spatial coordinates of the dataset elements
  # TODO: let the user determine the degree/order of the B-splines
  # the functions build.*.model need to return a list with 
  # an element B
  if(!is.null(spatial)) {
    
    if( is.null(spatial$autofill) ) {
      spatial$autofill = TRUE
    }
    
    # Splines model from Cappa & Cantet (2007)
    if(spatial$model == 'Cappa07') {
      sp <- build.splines.model(spatial$coord,
                                spatial$n.knots,
                                spatial$autofill,
                                degree = 3)
      effect.item <- list(name   = spatial$model,
                          pos    = pos - 1 + 1:ncol(sp$B),
                          levels = 1,
                          type   = 'cov',
                          model  = 'user_file_i',
                          file   = 'spatial',
                          var    = spatial$var.ini,
                          sp     = sp)
      
    }

    # Kronecker product of Autoregressive models
    # on the rows and columns (regular grids only)
    if(spatial$model == 'AR') {
      sp <- build.ar.model(spatial$coord,
                           spatial$rho,
                           spatial$autofill)
      # The number of levels of the effect must be the size
      # of the covariance matrix U. The Incidence matrix
      # need not contain an observation of the last level.
      stopifnot(identical(max(sp$U[, 1]), max(sp$U[, 2])))
      effect.item <- list(name   = spatial$model,
                          pos    = pos,
                          levels = max(sp$U[,1]),
                          type   = 'cross',
                          model  = 'user_file',
                          file   = 'spatial',
                          var    = spatial$var.ini,
                          sp     = sp)
      
    }
    
    # Blocks effect
    # An independent random effect
    # re-use the already detected random effect from the model-frame
    # change its name, and complete fields
    if(spatial$model == 'blocks') {
      sp <- build.blocks.model(spatial$coord,
                               as.numeric(spatial$id),
                               spatial$autofill)

      effect.item <- list(name   = spatial$model,
                          pos    = pos,
                          levels = nlevels(spatial$id),
                          type   = 'cross',
                          model  = 'diagonal',
                          file   = '',
                          var    = spatial$var.ini,
                          sp     = sp)
    }
    
    effects <- c(effects, 
                 spatial = list(effect.item))
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
    if(!is.null(fl)) {
    # NAs are written as 0
    write.table(fl$file, file=file.path(dir, fl$fname),
                row.names = FALSE, col.names = FALSE, na = "0")   
    # file.show(file.path(dir, fl$fname))
    }
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
  isF <- sapply(attr(mt, 'dataClasses'), 
                function(x) x %in% c('factor', 'ordered'))
  
  # Flags for specific effects
  isGenetic <- exists('genetic', as.environment(effects))
  isSpatial <- exists('spatial', as.environment(effects))
  
  # write labels for factor levels in results
  for( x in names(isF)[which(isF)] )
    rownames(result[[x]]) <- levels(mf[[x]])
  
  # Random and Fixed effects indices with respect to the 'effects' list
  fixed.effects.idx <- sapply(effects, function(x) !exists('model', x))
  diagonal.effects.idx <- sapply(effects, function(x) identical(x$model, 'diagonal'))
  special.effects.idx <- !(fixed.effects.idx | diagonal.effects.idx)
  random.effects.idx <- diagonal.effects.idx | special.effects.idx
  #   special.effects.idx <- which(names(effects) == 'genetic' | 
  #                                   names(effects)=='spatial')
  #   random.effects.idx <- c(diagonal.effects.idx, special.effects.idx)
  #   if( length(random.effects.idx) )
  #     fixed.effects.idx <- (1:length(effects))[-random.effects.idx]  else
  #       fixed.effects.idx <- (1:length(effects))
  #   diagonal.effects.idx <- which(attr(mt, 'term.types') == 'random')
  
  # Fixed effects coefficients
  beta <- sol.file$value[sol.file$effect %in% which(fixed.effects.idx)]
  #   .getXlevels(mt, mf)
  
  # Coefficients for model frame
  mf_values <- sol.file$value[sol.file$effect <= 
                                length(attr(mt, 'term.labels'))]
  
  # Random effects coefficients
  # TODO: Return Standard Errors as well.
  # How to compute standard errors of splines predicted values?
  ranef <- list()
  if( sum(diagonal.effects.idx) )
    ranef <- c(ranef, result[diagonal.effects.idx])
  
  genetic.fit <- spatial.fit <- genetic.pred <- spatial.pred <- NULL
  
  if(isGenetic){
    ranef$genetic <- result$genetic$value
    genetic.fit <- result$genetic$value[effects$genetic$idx]
  }
  # Spatial Surface
  # Here, the random effects are the underlying model parameters,
  # the fit is the predicted values for all the observations (which may include
  # several in the same location), and the pred is the predicted value in a full
  # rectangular grid, even if there were no observations there.
  if (isSpatial) {
    ranef$spatial <- result$spatial$value
    
    if( length(effects$spatial$pos) > 1 ){
      # Splines model
      spatial.fit <- data.frame(effects$spatial$sp$coord,
                                z = as.vector(effects$spatial$sp$B
                                              %*% 
                                                result$spatial$value))
      spatial.pred <- data.frame(effects$spatial$sp$plotting$grid,
                                 z = as.vector(effects$spatial$sp$plotting$B
                                               %*% result$spatial$value))
    } else if( effects$spatial$name == "AR" ) {
      # Autoregressive model
      # In the ordering of the dataset
      spatial.fit <- data.frame(effects$spatial$sp$coord,
                                z = result$spatial$value[effects$spatial$sp$B])
      spatial.pred <- data.frame(effects$spatial$sp$plotting$grid,
                                 z = result$spatial$value)
    } else if( effects$spatial$name == "blocks" ) {
      # Blocks model
      spatial.fit <- data.frame(effects$spatial$sp$coord,
                                z = result$spatial$value[effects$spatial$sp$B])
      tmp <- rep(NA, nrow(effects$spatial$sp$plotting$grid))
      tmp[effects$spatial$sp$map] <- result$spatial$value[effects$spatial$sp$B]
      spatial.pred <- data.frame(effects$spatial$sp$plotting$grid,
                                 z = tmp)
    }
  }
  # Build up the model matrix *for the fixed and random terms*
  # with one dummy variable per level of factors
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
  
  # This includes fixed and unstructured random effects
  eta <- drop(mm %*% mf_values)
  
  #   # Extract the BLUP for a given random effect result
  #   # unstructured random effects return a list with value and se
  #   # while spatial and genetic return a vector of values
  #   get_ranvalues <- function(x) {
  #     if(is.list(x)) return(x$value)
  #     else return(x)
  #   }
  #   gen.or.sp.idx <- names(ranef) %in% c('genetic', 'spatial')
  #   ranvalues <- lapply(ranef[gen.or.sp.idx], get_ranvalues)
  #   if(length(ranvalues)) 
  #   eta <- eta + rowSums(do.call(cbind, 
  #                                ranef[c('genetic', 'spatial')]))
  if(isGenetic | isSpatial)
    eta <- eta + rowSums(cbind(genetic.fit, spatial.fit$z))
  # Fitted Values
  # ASSUMPTION: Linear Model (not generalized)
  # TODO: apply inverse link
  mu = eta

  # REML info
  # TODO: 'delta convergence' only for AI-REML?
  reml.ver <- sub('^\\s+([[:graph:]]* +ver\\. +[0-9.]*).*$', '\\1', 
                  grep('REML', reml.out, value = TRUE))
  last.round.idx <- tail(grep('In round', reml.out), 1)
  last.round <- as.numeric(strsplit(strsplit(reml.out[last.round.idx],
                                             split='In round')[[1]][2],
                                    split='convergence=')[[1]])
  
  # Maximum number of iterations
  # (Hardcoded in REML and AIREML)
  max.it <- 5000
  
  # Variance components
  if( identical(last.round[1], max.it) ) {
    warning('The algorithm did not converge')
    varcomp <- cbind('Estimated variances' = rep(NA, sum(random.effects.idx) + 1L))
    rownames(varcomp) <- c(names(effects)[random.effects.idx], 'Residual')
  } else {
    # Variance components
    # ASSUMPTION: I always have a Genetic Variance component
    
    # Issue #2 In Linux, AIREMLF90 prints S.D. for R and G
    # while under Windows it outputs SE for R and G
    # Update: from version 1.109 (at least), Linux updated to SE as well
    #   sd.label <- ifelse(.Platform$OS.type == 'windows', 'SE', 'S.D.')
    sd.label <- ifelse(TRUE, 'SE', 'S.D.')
    
    varcomp.idx <- grep('Genetic variance|Residual variance', reml.out) + 1
    # There should be one variance for each random effect plus one resid. var.
    stopifnot(identical(length(varcomp.idx), sum(random.effects.idx) + 1L))
    varcomp <- cbind('Estimated variances' = as.numeric(reml.out[varcomp.idx]))
    rownames(varcomp) <- c(names(effects)[random.effects.idx], 'Residual')
    
    # EM-REML does not print Standard Errors for variance components
    if(method == 'ai'){
      varsd.idx <- grep(paste(sd.label, 'for G|for R'), reml.out) + 1
      # There should be one variance for each random effect plus one resid. var.
      stopifnot(identical(length(varcomp.idx), sum(random.effects.idx) + 1L))
      
      # Watch out!! AI-REML gives SE for R in the *first place*
      # even when the variance component was last
      varsd.idx <- c(varsd.idx[-1], varsd.idx[1])
      varcomp <- cbind(varcomp, 'S.E.' = as.numeric(reml.out[varsd.idx]))
    }

    # Update: we step back from this decision. We report the estimated variance
    # parameter of the spatial effect, even if it is not additive.
    # It might be used for comparison of the same model fitted to different data.
    #     # If spatial, report the observed spatial variance, rather than the
    #     # estimated variance of the spline effects which is meaningless
    #     # TODO: Return the true field variance multiplying matrices and whatever
    #     if(isSpatial) { 
    #       varcomp['spatial', 1] <- var(ranef$spatial)
    #       if(method == 'ai') varcomp['spatial', 2] <- NA
    #     }
    
    # For additive variance decomposition, we return as well the covariance matrices
    # of each component.
    # V(y) = \sigma_u^2 Z A Z' + \sigma_v^2 B U B' + \sigma_e^2 I
    # This is not giving me anything additive ???
  }
  
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
    components = list(pedigree = isGenetic,
                      spatial  = isSpatial), # TODO competition, ...
    effects = effects,
    mf = mf,
    mm = mm,
    y = y,
    fixed = result[fixed.effects.idx],
    ranef = ranef,
    eta = eta,
    mu = mu,
    residuals = y - mu,
    genetic = list(fit        = genetic.fit),
    spatial = list(name       = effects$spatial$name,
                   model      = effects$spatial$sp[1:3],
                   fit        = spatial.fit,
                   prediction = spatial.pred),
    var = varcomp,
    fit = fit,
    reml = reml
  )
  return(ans)
}


#' Build Model Frame
#' 
#' Merges fixed and random terms into a single call
#' and returns the corresponding model frame
#' optionally removing the intercept term
build.mf <- function(call) {
	terms.list <- list()
	
  # Don't use the model.frame intercept
  # Make sure there is something for the rhs
  terms.list$int <- 0

  ## Fixed effects
  fxd <- eval(call$fixed, parent.frame(2))
  tfxd <- terms(fxd)
  
	# Add an intercept manually only if the user requested it
  # *and* there are no other *categorical* fixed effects
  tempmf <- eval(call('model.frame',
                      formula = fxd,
                      data = quote(data)),
                 parent.frame())
  tempt <- terms(tempmf)
  tempc <- attr(tempt, 'dataClasses')[attr(tempt, 'term.labels')]
  any.cat <- any(tempc %in% c('factor', 'ordered'))
	
  if( attr(tfxd, 'intercept') == 1L & !any.cat ) {
	  fxd <- update(fxd, " ~ Intercept + .")
	}
	terms.list$fxd <- attr(terms(fxd), 'term.labels')
	
	
	## Random effects (unstructured)
	rnd <- eval(call$random, parent.frame(2))
  if(!is.null(rnd))
    terms.list$rnd <- attr(terms(rnd), 'term.labels')
	
	## Join fixed and random
	lhs <- as.character(fxd[[2]])
	rhs <- paste(do.call(c, terms.list), collapse = '+')  
	fml <- as.formula(paste(lhs, rhs, sep = '~'), env = parent.frame(2))
	
  # Build Model Frame
  # Use na.pass to allow missing observations which will be handled later
	mfcall <- call('model.frame',
                 formula = fml,
                 data = quote(transform(data, Intercept = 1)),
                 na.action = na.pass)
	mf <- eval(mfcall, parent.frame())
  mt <- attr(mf, 'terms')
  
  
  # Add attribute indicating 'fixed' or 'random'
  stopifnot(length(attr(mt, 'term.labels')) == 
              length(terms.list$fxd) + length(terms.list$rnd))
  label_var <- function(x, label) sapply(x, function(x) label)
  tl <- c(label_var(terms.list$fxd, 'fixed'),
          label_var(terms.list$rnd, 'random'))
  attr(attr(mf, 'terms'), 'term.types') <- tl

	
	## Strings as factors
	str.idx <- which(attr(mt, 'dataClasses') == 'character' | 
	                   attr(mt, 'dataClasses') == 'other')
	if(length(str.idx)) {
	  mf[str.idx] <- lapply(mf[str.idx], as.factor)
    attr(attr(mf, 'terms'), 'dataClasses')[str.idx] <- 'factor'
	}
  return(mf)
}