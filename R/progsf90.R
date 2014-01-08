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

