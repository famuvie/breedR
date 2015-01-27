# progsf90 class
# 
# This function parses a model frame and extracts the relevant fields
# that are to be written in the parameter, data and auxiliary files
# of the progsf90 programs.
progsf90 <- function (mf, effects, opt = c("sol se"), res.var.ini = 10) {
  
  
  # Build models for random effects
  mt <- attr(mf, 'terms')
  random.effects.idx <- c(which(attr(mt, 'term.types') == 'random'),
                          which(names(effects) %in% c('genetic', 'pec', 'spatial')))
  
                          
  parse.rangroup <- function(x) {
    group.size <- nrow(as.matrix(effects[[x]]$var))
    group.head <- head(which(effects[[x]]$levels != 0), 1)
    # Determine the right position in the effects list
    group.head.abs <- sum(sapply(effect.lst, length)[1:(x-1)]) + group.head
    return(list(pos = group.head.abs + 1:group.size - 1,
                type = effects[[x]]$model, 
                file = effects[[x]]$file, 
                cov  = effects[[x]]$var))
  }
  # Number of traits
  # (size of the response vector or matrix)
  ntraits <- ncol(as.matrix(model.response(mf)))
  
  # Weights
  weights <- ''     # No weights for the moment --- TODO
  
  # Builds the lines in the EFFECTS section
  effect.lst <- sapply(effects,
                       function(x) with(x, paste(pos, levels, type)))
  # Parameters  
  par <- list(datafile = 'data',
              ntraits  = ntraits,
              neffects = sum(sapply(effects, function(x) length(x$pos))),
              observations = 1:ntraits,
              weights  = weights,
              effects  = effect.lst,
              residvar = res.var.ini,
              rangroup = lapply(random.effects.idx, parse.rangroup),
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
           genetic       = effects$genetic$gen$B,
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


# Build effects parameters
# 
# This function builds a list of effects parameters
# as required by Misztal's progsf90 suite of programs
# @references
# \url{http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90.pdf}
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
  if( attr(mt, 'intercept') ) {
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
  
  # Parameters for all effects in the fixed formula
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
  # Both the additive genetic and spatial terms are "cross"
  # The competition effect is nested into the additive genetic and it is 'cov'
  # with the funny structure of zeroes to produce the sum of effects
  # For pedigree effects, there might be more levels than those
  # present in the data. We should declare the levels present in the pedigree.
  
  if( !is.null(genetic) ) {

    if( is.null(genetic$autofill) ) {
      genetic$autofill = TRUE
    }
    
    gen <- build.genetic.model(genetic)
    
    # Initial values for the add-animal model, which
    # are present also in the competition model
    gen.levels <- nrow(as.data.frame(genetic$pedigree))
    gen.type   <- 'cross'

    # max n of competitors (0 under no competition)
    n.comp <- (ncol(gen$B) - 1)/2
    
    # In the competition model, B has one column for the animal direct effect,
    # and two columns for each neighbouring effect (at most 8), one for the 
    # IC and the other for keeping the neigbour index.
    if( genetic$model == 'competition' ) {
      stopifnot(n.comp > 0)
      gen.levels <- c(gen.levels,
                      rep(0, n.comp - 1),
                      gen.levels)
      gen.type   <- c(gen.type,
                      paste('cov', pos + n.comp + 1:n.comp))
    }

    effects <- c(effects, 
                 genetic = list(
                   list(pos    = pos - 1 + 1:(1 + n.comp),
                        levels = gen.levels,
                        type   = gen.type,
                        model  = 'add_animal', # both add_animal or competition
                        file   = 'pedigree',
                        ped    = as.data.frame(genetic$pedigree),
                        var    = genetic$var.ini,
                        gen    = gen)))
    
    # Permanent Environmental Competition Effect
    # In this case, create yet another effect, with the same incidence matrix
    # but unstructured
    if( genetic$model == 'competition' ) {
      if( genetic$pec$present )
        effects <- c(effects, 
                     pec = list(
                       list(pos    = pos + 1:n.comp,
                            levels = gen.levels[-1],
                            type   = gen.type[-1],
                            model  = 'diagonal', # both add_animal or competition
                            file   = '',
                            var    = genetic$pec$var.ini)))
    }
    
    pos = pos + 1 + 2 * n.comp   # under add_animal, n.comp = 0
  }
  
  # Spatial effect if applicable
  # We only have spatial coordinates of the dataset elements
  # TODO: let the user determine the degree/order of the B-splines
  # the functions build.*.model need to return a list with 
  # an element B
  if( !is.null(spatial) ) {
    
    if( is.null(spatial$autofill) ) {
      spatial$autofill = TRUE
    }
    
    # Splines model from Cappa & Cantet (2007)
    if(spatial$model == 'splines') {
      sp <- build.splines.model(spatial$coord,
                                spatial$n.knots,
                                spatial$autofill,
                                degree = 3)
      n.pos <- ncol(sp$B)
      effect.item <- list(name   = spatial$model,
                          pos    = pos - 1 + 1:n.pos,
                          levels = c(rep(0, n.pos-1), n.pos),
                          type   = paste('cov', pos - 1 + n.pos + 1:n.pos),
                          model  = 'user_file_i',
                          file   = 'spatial',
                          var    = spatial$var.ini,
                          sp     = sp)
      pos = pos + 2*ncol(sp$B)
      
    }

    # Kronecker product of Autoregressive models
    # on the rows and columns (regular grids only)
    if(spatial$model == 'AR') {
      sp <- build.ar.model(spatial$coord,
                           spatial$rho,
                           spatial$autofill)
      # The number of levels of the effect must be the size
      # of the covariance/precision matrix U. The Incidence matrix
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
      pos = pos + 1
    }
    
    # Blocks effect
    # An independent random effect
    # re-use the already detected random effect from the model-frame
    # change its name, and complete fields
    if(spatial$model == 'blocks') {
      # spatial$id is always a factor at this point
      stopifnot(is.factor(spatial$id))
      
      sp <- build.blocks.model(spatial$coord,
                               spatial$id,
                               spatial$autofill)

      effect.item <- list(name   = spatial$model,
                          pos    = pos,
                          levels = nlevels(spatial$id),
                          type   = 'cross',
                          model  = 'diagonal',
                          file   = '',
                          var    = spatial$var.ini,
                          sp     = sp)
      pos = pos + 1
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
           c(sapply(rangroup, 
                  function(x) c('RANDOM_GROUP',
                                paste(x$pos, collapse = ' '),
                                'RANDOM_TYPE',
                                x$type,
                                'FILE',
                                x$file,
                                '(CO)VARIANCES',
                                apply(as.matrix(x$cov), 1,
                                      paste,
                                      collapse = ' '))),
             recursive = TRUE),
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



# Parse results from a progsf90 'solutions' file
parse_results <- function (solfile, effects, mf, reml.out, method, mcout) {
  
  ## Parsing the results
  sol.file <- read.table(solfile, header=FALSE, skip=1)
  colnames(sol.file) <- c('trait', 'effect', 'level', 'value', 's.e.')
  
  # Assuming one trait only
  result <- by(sol.file[,4:5], sol.file$effect, identity)
  
  # Different results can be associated to a single (group) effect
  effect.size <- function(x) ifelse(exists('var', x),
                                   nrow(as.matrix(x$var)), 1)
  
  # for each result, a pointer to its corresponding effect
  result_effect.map <- 
    unlist(sapply(seq_along(effects),
                  function(idx) rep(idx, effect.size(effects[[idx]]))))
    
  # Name the results according to effects
  # Effects can be grouped (e.g. competition) and account for correlated
  # effects
  names(result) <- names(effects)[result_effect.map]
  
  # Account for competition models
  # TODO: Do this more generally, as more effects can be grouped
  if ( sum(names(result) == 'genetic') == 2 ) {
    names(result)[names(result) == 'genetic'] <- c('genetic-direct',
                                                   'genetic-competition')
  }
  
  
  # Identify factors in model terms
  mt <- attr(mf, 'terms')
  isF <- sapply(attr(mt, 'dataClasses'), 
                function(x) x %in% c('factor', 'ordered'))
  
  # Flags for specific effects
  isGenetic <- exists('genetic', effects)
  isSpatial <- exists('spatial', effects)
  
  # write labels for factor levels in results
  for( x in names(isF)[which(isF)] )
    rownames(result[[x]]) <- levels(mf[[x]])
  
  # Random and Fixed effects indices with respect to the 'effects' list
  fixed.effects.idx <- sapply(effects, function(x) !exists('model', x))
  diagonal.effects.idx <- sapply(effects,
                                 function(x) identical(x$model, 'diagonal'))
  special.effects.idx <- !(fixed.effects.idx | diagonal.effects.idx)
  random.effects.idx <- diagonal.effects.idx | special.effects.idx
  

  
  # Random effects coefficients
  # TODO: Return Standard Errors as well.
  # How to compute standard errors of splines predicted values?
  ranef <- result[result_effect.map %in% which(random.effects.idx)]
  #   ranef <- list()
  #   if( sum(diagonal.effects.idx) )
  #     ranef <- c(ranef, result[diagonal.effects.idx])
  



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
    sd.label <- ifelse(TRUE, 'SE', 'S.D.')
    
    varcomp.idx <- grep('Genetic variance|Residual variance', reml.out) + 1
    # There should be one variance for each random effect plus one resid. var.
    stopifnot(identical(length(varcomp.idx), sum(random.effects.idx) + 1L))

    rangroup.sizes <- c(sapply(effects[random.effects.idx],
                               effect.size),
                        resid = 1)

    if( all(rangroup.sizes == 1) ){
      varcomp <- cbind('Estimated variances' = as.numeric(reml.out[varcomp.idx]))
      rownames(varcomp) <- c(names(effects)[random.effects.idx], 'Residual')
    } else {
      varcomp.str <- lapply(mapply(function(x, y) x + 1:y,
                                   varcomp.idx-1,
                                   rangroup.sizes),
                            function(x) reml.out[x])
      parse.varmat <- function(v, names) {
        # get the numeric values from the strings, spliting by spaces
        ans <- sapply(v, function(x) as.numeric(strsplit(x, ' +')[[1]][-1]))
        # if no sub-names passed, remove all naming
        if( is.null(names) ) names(ans) <- dimnames(ans) <- NULL
        # otherwise, put the given names in both dimensions (covariance matrix)
        else dimnames(ans) <- list(names, names)
        ans
      }

      # names for the members of a group (if more than one)
      get_subnames <- function(name) {
        ran.idx <- grep(name, names(ranef))
        if( length(ran.idx) > 1 ) {
          sn <- sapply(names(ranef)[ran.idx],
                       function(x) strsplit(x, paste0(name, '-'))[[1]][2])
        } else sn <- NULL
        return(sn)
      }
      
      subnames <- sapply(names(rangroup.sizes), get_subnames)
      varcomp <- mapply(parse.varmat, varcomp.str, subnames)
      names(varcomp) <- c(names(effects)[random.effects.idx], 'Residual')
    }
    
    # EM-REML does not print Standard Errors for variance components
    if(method == 'ai'){
      varsd.idx <- grep(paste(sd.label, 'for G|for R'), reml.out) + 1
      # There should be one variance for each random effect plus one resid. var.
      stopifnot(identical(length(varcomp.idx), sum(random.effects.idx) + 1L))
      
      # Watch out!! AI-REML gives SE for R in the *first place*
      # even when the variance component was last
      varsd.idx <- c(varsd.idx[-1], varsd.idx[1])
      if( all(rangroup.sizes == 1) ){
        varcomp <- cbind(varcomp, 'S.E.' = as.numeric(reml.out[varsd.idx]))
      } else {
        varsd.str <- lapply(mapply(function(x, y) x + 1:y,
                                     varsd.idx-1,
                                     rangroup.sizes),
                              function(x) reml.out[x])
        varsd <- mapply(parse.varmat, varsd.str, subnames)
        names(varsd) <- c(names(effects)[random.effects.idx], 'Residual')
      }
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
    fixed = result[which(fixed.effects.idx)],
    ranef = ranef,
#     eta = eta,
#     mu = mu,
#     residuals = y - mu,
#     genetic = list(fit        = genetic.fit),
#     spatial = list(name       = effects$spatial$name,
#                    model      = effects$spatial$sp[1:3],
#                    fit        = spatial.fit,
#                    prediction = spatial.pred),
    var = varcomp,
    fit = fit,
    reml = reml
  )
  return(ans)
}


# Build Model Frame
# 
# Merges fixed and random terms into a single call
# and returns the corresponding model frame
# optionally removing the intercept term
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