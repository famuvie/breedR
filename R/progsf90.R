# Build effects parameters
# 
# This function builds a list of effects parameters
# as required by Misztal's progsf90 suite of programs
# @references
# \url{http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90.pdf}
build.effects <- function (mf, genetic, spatial, generic, var.ini) {
  
  # Build up effects data (position, levels, type)
  
  # Model terms
  mt <- attr(mf, 'terms')
  
  #   # Number of traits
  #   # (size of the response vector or matrix)
  #   ntraits <- ncol(as.matrix(stats::model.response(mf)))
  
  #   # Position counter
  #   pos = ntraits + 1
  
  effects <- list()
  
  #   # Intercept term, if not precluded explicitly in the formula
  # This is taken care of in build.mf()
  #   if (attr(mt, 'intercept')) {
  #     effects <- c(effects, list('(Intercept)' = fixed(rep(1, nrow(mf)))))
  # #     effects <- list('(Intercept)'=list(pos = pos, levels=1L, type='cross'))
  # #     pos <- pos + 1
  #   }
  
  ## Fixed effects
  effects <- c(effects,
               lapply(mf[, names(which(attr(mt, 'term.types') == 'fixed')),
                         drop = FALSE],
                      fixed))
  
  ## Random effects
  mf.rnd <- mf[, names(which(attr(mt, 'term.types') == 'random')),
               drop = FALSE]

  for (x in names(mf.rnd)) {
    effect.item <- effect_group(list(diagonal(mf.rnd[[x]])),
                                cov.ini = var.ini[[x]],
                                ntraits = ncol(stats::model.response(mf)))
    effect.item.list <- structure(list(effect.item),
                                  names = x)
    effects <- c(effects, effect.item.list)
  }
  
  # Genetic effect
  # Both the additive genetic and spatial terms are "cross"
  # The competition effect is nested into the additive genetic and it is 'cov'
  # with the funny structure of zeroes to produce the sum of effects
  # For pedigree effects, there might be more levels than those
  # present in the data. We must declare the levels present in the pedigree.
  
  if (!is.null(genetic)) {
    
    ## additive-genetic effect
    gen_direct <- additive_genetic_animal(
      pedigree = genetic$pedigree,
      idx      = genetic$id)
    
    ## effect-group with direct and eventually competition effects
    if (genetic$model == 'add_animal') {
      
      ## additive-genetic only
      effect.item <- effect_group(list(direct = gen_direct),
                                  cov.ini = genetic$var.ini,
                                  ntraits = ncol(stats::model.response(mf)))
      effect.item.list <- list(genetic = effect.item)
    } else {
      
      ## both direct and competition effects
      stopifnot(genetic$model == 'competition')
      gen_comp <- additive_genetic_competition(
        pedigree    = genetic$pedigree,
        coordinates = genetic$coordinates,
        id          = genetic$id,
        decay       = genetic$competition_decay,
        autofill    = genetic$autofill)
      
      effect.item <- effect_group(list(genetic_direct = gen_direct,
                                       genetic_competition = gen_comp),
                                  cov.ini = genetic$var.ini,
                                  ntraits = ncol(stats::model.response(mf)))
      effect.item.list <- list(genetic = effect.item)
      
      ## eventually, a second effect-group for pec      
      if (genetic$pec$present) {
        pec <- permanent_environmental_competition(
          coordinates = genetic$coordinates,
          decay       = genetic$competition_decay,
          autofill    = genetic$autofill
        )
        effect.item <- effect_group(list(pec = pec),
                                    cov.ini = genetic$pec$var.ini,
                                    ntraits = ncol(stats::model.response(mf)))
        effect.item.list <- c(
          effect.item.list,
          list(pec = effect.item)
        )
      }
    }
    
    effects <- c(effects, effect.item.list)
  }
  
  # Spatial effect if applicable
  # We only have spatial coordinates of the dataset elements
  # TODO: let the user determine the degree/order of the B-splines
  # the functions build.*.model need to return a list with 
  # an element B
  if (!is.null(spatial)) {
    
    ## call the corresponding model constructor
    breedr_builder <- paste('breedr', tolower(spatial$model), sep = '_')
    sp <- do.call(breedr_builder, spatial)
    
    ## build the effect group with the (single) spatial model
    effect.item <- effect_group(structure(list(sp), names = class(sp)[1]),
                                cov.ini = spatial$var.ini,
                                ntraits = ncol(stats::model.response(mf)))
    
    effects <- c(effects, list(spatial = effect.item))
  }
  
  
  ## Generic effect if applicable
  
  if( !is.null(generic) ) {
    
    ## From each element in the generic list, we build a generic object
    ## and make an effect_group with it alone, and the corresponding var.ini
    make_group <- function(x) {
      stopifnot('var.ini' %in% names(x))
      go <- do.call('generic', x[-grep('var.ini', names(x))])
      ef <- effect_group(list(go),
                         x[['var.ini']],
                         ntraits = ncol(stats::model.response(mf)))
      return(ef)
    }
    generic.groups <- lapply(generic, make_group)
    
    ## Make sure the names do not clash any of the 'special' names
    ## i.e. 'genetic' or 'spatial'
    match.idx <- names(generic.groups) %in% c('genetic', 'spatial')
    if (any(match.idx)) {
      names(generic.groups)[match.idx] <- 
        paste0('generic_', names(generic.groups)[match.idx])
    }
    
    effects <- c(effects,
                 generic.groups)
  }
  
  return(effects)
}


#' progsf90 class
#' 
#' This function parses a model frame and extracts the relevant fields
#' that are to be written in the parameter, data and auxiliary files
#' of the progsf90 programs.
#' 
#' @param mf model.frame for fixed and diagonal random effects
#' @param weights a vector of weights for the residual variance
#' @param effects breedr_modelframe
#' @param opt character. Options to be passed to Misztal's programs
#' @param res.var.ini positive number. Initial value for the residual variance.
#' @family progsf90
progsf90 <- function (mf, weights, effects, opt = c("sol se"), res.var.ini = 10) {
  
  ## Build models for random effects (in 'random')
  mt <- attr(mf, 'terms')
  
  ## types of effects (fixed / random)
  eff_types <- sapply(effects, effect_type)
  
  ## names of random groups
  rangroup.idx <- which(eff_types == 'random')
  
#   random.effects.idx <- 
#     unique(c(which(attr(mt, 'term.types') == 'random'),
#       ## This only works temporarily. After completing the refactoring,
#       ## I should use the method effect_type() for each element in effect.
#       which(sapply(effects, inherits, 'effect_group'))))
  
  # Number of traits
  # TODO: Should I pass only the response instead of the whole mf?
  # or maybe only ntraits?, or maybe include the response in effects?
  # (size of the response vector or matrix)
  ntraits <- ncol(as.matrix(stats::model.response(mf)))
  
  # Weights position
  w_pos <- ifelse(is.null(weights), '', ntraits + 1)

  ## renderpf90 all the effects
  ## positions in data file starting after offset for traits and weights
  effects.pf90 <- renderpf90.breedr_modelframe(effects, ntraits, w_pos>0)

  # Builds the lines in the EFFECTS section
  
  setup_effectline <- function(x) {
    with(x, paste(pos, levels, type, nest))
  }
  
  effect.lst <- lapply(effects.pf90, setup_effectline)

  # Phenotype
  Y <- mf[, attr(attr(mf, 'terms'), 'response')]
  
  # Code for missing observations
  # Use 0 if outside range of variation. Otherwise, use option 'missing'.
  # Overriden by user specification
  if (!any(grepl('missing', opt))) {
    missing_code <- pf90_code_missing(Y)
    Y[which(is.na(Y))] <- missing_code
    if (!identical(missing_code, 0)) {
      opt <- c(opt, paste('missing', missing_code))
    }
  }
  
  parse.rangroup <- function(x) {
    group.size <- nrow(as.matrix(effects.pf90[[x]]$var))/ntraits
    ## The group 'head' is the first effect with a number of levels > 0
    group.head <- head(which(effects.pf90[[x]]$levels != 0), 1)
    # Determine the right position in the effects list
    group.head.abs <- sum(sapply(effect.lst, length)[1:(x-1)]) + group.head
    
    return(list(pos = group.head.abs + 1:group.size - 1,
                type = effects.pf90[[x]]$model, 
                file = effects.pf90[[x]]$file_name, 
                cov  = effects.pf90[[x]]$var))
  }
  # Parameters  
  par <- list(datafile = 'data',
              ntraits  = ntraits,
              neffects = sum(sapply(effect.lst, length)),
              observations = paste(seq_len(ntraits), collapse = " "),
              weights  = w_pos,
              effects  = effect.lst,
              residvar = res.var.ini,
              rangroup = lapply(rangroup.idx, parse.rangroup),
              options = opt
  )
  
  # Data 
  # Columns ordered as in the effects list
  # after the trait(s) 
  dat.l <- lapply(effects.pf90, function(x) x$data)
  if (length(unique(nr <- sapply(dat.l, nrow))) > 1) {
    most.common.size <- as.numeric(names(which.max(table(nr))))
    offending.effects <- which(nr != most.common.size)
    msg <- paste('The incidence matrix of',
                  toString(names(offending.effects)),
                  'has not the expected number of rows.')
    stop(msg)
  }
  dat <- do.call(cbind, 
                 c(list(phenotype = Y),
                   list(w = weights),
                   dat.l))
  
  # Forbid missing values in dependent variables
  if (anyNA(dat[, -(1:ntraits)])) {
    idx <- which(is.na(dat[, -(1:ntraits), drop = FALSE]), arr.ind = TRUE)[, 1]
    stop('\nMissing values in covariates are not allowed\n',
         paste('check individuals:', paste(idx, collapse = ', ')))
  }

  # Additional Files
  build.file.single <- function(ef) {
    
    if (ef$file_name != '') {
      return(list(fname = ef$file_name,
                  file  = ef$file))
    }
  }
  files <- lapply(effects.pf90[rangroup.idx], build.file.single)
  
  ans <- list(parameter = par, data = dat, files = files)
  class(ans) <- 'progsf90'
  return(ans)
}


write.progsf90 <- function (pf90, dir) {
  
  write_matrix <- function(x) {
    apply(as.matrix(x), 1, paste, collapse = " ")
  }
  
  write_rangroup <- function(x) {
    c('RANDOM_GROUP',
      paste(x$pos, collapse = ' '),
      'RANDOM_TYPE',
      x$type,
      'FILE',
      x$file,
      '(CO)VARIANCES',
      write_matrix(x$cov))
  }
  
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
           'RANDOM_RESIDUAL VALUES', write_matrix(residvar),
           c(sapply(rangroup, write_rangroup), recursive = TRUE),
           paste('OPTION', options)))
  
  writeLines(as.character(parameter.file),
             con = file.path(dir, 'parameters'))
  # file.show(parameter.file.path)
  
  # Data file
  utils::write.table(pf90$data,
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

  
  ## Parse the AI matrix from AIREML output
  parse_invAI <- function(x) {
    ## precedent and subsequent line numbers
    idx <- grep('inverse of AI matrix', x)
    stopifnot(length(idx) == 2)
    mat.txt <- x[(idx[1]+1):(idx[2]-1)]
    
    invAI <- parse.txtmat(mat.txt)
    return(invAI)
  }

  
  ## Parse functions of (co)variances
  ## from OPTION se_covar_function
  parse_functions <- function(x) {
    
    parse_function <- function(x) {
      pick_nmbr <- function(x) {
        conv <- suppressWarnings(as.numeric(x))
        stopifnot(sum(idx <- !is.na(conv)) == 1)
        return(conv[idx])
      }
      structure(sapply(strsplit(x, ' '), pick_nmbr),
                names = c('mean', 'sample mean', 'sample sd'))
    }
    
    pattern <- "^.*? SE for function of \\(co\\)variances (\\S+) .*$"
    labels.idx <- grep(pattern, x)
    labels <- gsub(pattern, "\\1", x[labels.idx])
    idx <- grep('  - Function:', x)
    stopifnot((N <- length(idx)) == length(labels.idx))
    
    ans <- sapply(lapply(idx, function(i) x[i+1:3]), parse_function)
    if (N > 0L) colnames(ans) <- labels
    return(ans)
  }

  
  ## Number and names of traits
  ntraits <- as.numeric(tail(unlist(strsplit(
    grep("Number of Traits", reml.out, value = TRUE),
    ' +')), 1))
  trait_names <- colnames(stats::model.response(mf))  # NULL for 1 trait
  
  # Parsing the results
  sol.file <- try(utils::read.table(solfile, header=FALSE, skip=1))
  if( inherits(sol.file, 'try-error') ) {
    ## The output file is formatted with fixed-width columns
    ## leaving 4 spaces between the columns trait and effect
    ## If there are more than 999 total random effects (e.g.
    ## in a splines model with 30x30 knots), this two columns
    ## become one.
    ## The following is a workaround treating it as one column
    ## filling the last with NA and then rearranging
    ## Fixes #24
    sol.file <- read.table(solfile, header=FALSE, skip=1, fill=TRUE)
    bug.idx <- which(is.na(sol.file[, 5]))
    
    ## check that only one column was lost
    stopifnot(all(!is.na(sol.file[bug.idx, -5])))
    
    ## check that the first digit of the trait/effect code is 1
    stopifnot(all(substr(sol.file[bug.idx, 1], 1, 1) == "1"))
    
    ## correct the code of the effect by removing the leading 1
    sol.file[bug.idx, 1] <- sol.file[bug.idx, 1] - 1e4
    
    ## shift all the columns to the right
    sol.file[bug.idx, ] <- cbind(1, sol.file[bug.idx, 1:4])
  }
  colnames(sol.file) <- c('trait', 'effect', 'level', 'value', 's.e.')
  
  # trait < level < effect
  split_by <- function(x, var) {
    col.id <- match(var, names(x))
    split(x[, -col.id], x[[col.id]])
  }
  result_by_effect <- split_by(sol.file[, -3], 'effect')
  result <- lapply(result_by_effect, split_by, 'trait')

  # Name the results according to effects
  # Effects can be grouped (e.g. competition) and account for correlated
  # effects
  names(result) <- get_efnames(effects)
  
  # Name traits within effects
  result <- lapply(result, structure, names = trait_names)

  # Different results can be associated to a single (group) effect
  effect.size <- vapply(effects, dim, numeric(2))["size", ]

  # count fixed effects as of size 1
  effect.size[effect.size == 0] <- 1
  
  # for each result, a pointer to its corresponding effect
  result_effect.map <- 
    unlist(sapply(seq_along(effects),
                  function(idx) rep(idx, effect.size[idx])))
    
  # Identify factors in model terms
  mt <- attr(mf, 'terms')
  isF <- sapply(attr(mt, 'dataClasses'), 
                function(x) x %in% c('factor', 'ordered'))
  
  # Flags for specific effects
  isGenetic <- exists('genetic', effects)
  isSpatial <- exists('spatial', effects)
  
  # write labels for factor levels in results
  for (x in names(isF)[which(isF)] ) {
    for (trait in seq_len(ntraits)) {
      rownames(result[[x]][[trait]]) <- levels(mf[[x]])
    }
  }
  
  # Random and Fixed effects indices with respect to the 'effects' list
  effect.type <- vapply(effects, effect_type, '')
  

  # Random effects coefficients
  ranef <- result[result_effect.map %in% which(effect.type == 'random')]

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

  ## Dimension of random effects
  rangroup.sizes <- c(effect.size[effect.type == 'random'],
                      resid = 1)*ntraits

  ## index of variance component results
  varcomp.idx <- grep('Genetic variance|Residual variance', reml.out) + 1
  
  # Variance components
  if (identical(last.round[1], max.it)) {
    warning('The algorithm did not converge')
    varcomp <- cbind('Estimated variances' = rep(NA, sum(effect.type == 'random') + 1L))
    rownames(varcomp) <- c(names(effects)[effect.type == 'random'], 'Residual')
  } else {
    # Variance components
    sd.label <- ifelse(TRUE, 'SE', 'S.D.')
    
    # There should be one variance for each random effect plus one resid. var.
    stopifnot(identical(length(varcomp.idx), sum(effect.type == 'random') + 1L))

    if (all(rangroup.sizes == 1)){
      varcomp <- cbind('Estimated variances' = as.numeric(reml.out[varcomp.idx]))
      rownames(varcomp) <- c(names(effects)[effect.type == 'random'], 'Residual')
    } else {
      varcomp.str <- lapply(mapply(function(x, y) x + 1:y,
                                   varcomp.idx-1,
                                   rangroup.sizes,
                                   SIMPLIFY = FALSE),
                            function(x) reml.out[x])

      # names for the members of a group (if more than one)
      get_subnames <- function(name) {
        ## effect sub-names (or NULL)
        esn <- names(effects[[name]]$effects)
        names_effect(esn, trait_names)
      }
      
      subnames <- lapply(names(rangroup.sizes), get_subnames)
      varcomp <- mapply(parse.txtmat, varcomp.str, subnames, SIMPLIFY = FALSE)
      names(varcomp) <- c(names(effects)[effect.type == 'random'], 'Residual')
    }
    
    # EM-REML does not print Standard Errors for variance components
    if (method == 'ai') {
      varsd.idx <- grep(paste(sd.label, 'for G|for R'), reml.out) + 1
      # There should be one variance for each random effect plus one resid. var.
      stopifnot(identical(length(varcomp.idx), sum(effect.type == 'random') + 1L))
      
      if (all(rangroup.sizes == 1)){
        varcomp <- cbind(varcomp, 'S.E.' = as.numeric(reml.out[varsd.idx]))
      } else {
        varsd.str <- lapply(mapply(function(x, y) x + 1:y,
                                   varsd.idx-1,
                                   rangroup.sizes,
                                   SIMPLIFY = FALSE),
                            function(x) reml.out[x])
        varsd <- mapply(parse.txtmat, varsd.str, subnames, SIMPLIFY = FALSE)
        names(varsd) <- c(names(effects)[effect.type == 'random'], 'Residual')
        varcomp <- cbind("Estimated variances" = varcomp,
                         "S.E." = varsd)
      }
    }
  }
  
  # REML algorithm
  reml <- list(
    version = gsub('\\s+', ' ', reml.ver),
    rounds = last.round[1],
    convergence = last.round[2],
    delta.conv = as.numeric(strsplit(reml.out[last.round.idx+1],
                                     split='delta convergence=')[[1]][2]),
    output = reml.out
  )
  
  ## AI matrix
  if (method == 'ai') {

    comp_names <- unlist(
      mapply(vcnames, names(rangroup.sizes), rangroup.sizes,
             MoreArgs = list(trnames = trait_names), SIMPLIFY = FALSE))
    reml$invAI <- parse_invAI(reml.out)
    dimnames(reml$invAI) <- list(comp_names, comp_names)
  }

  # Fit info
  last.fit <- as.numeric(strsplit(strsplit(reml.out[last.round.idx-1],
                                           split='-2logL =')[[1]][2],
                                  split=': AIC =')[[1]])
  fit <- list(
    '-2logL' = last.fit[1],
    AIC = last.fit[2]
  )

  ans <- list(
    call = mcout,
    method = method,
    components = list(pedigree = isGenetic,
                      spatial  = isSpatial),
    effects = effects,
    mf = mf,
    fixed = result[result_effect.map %in% which(effect.type == 'fixed')],
    ranef = ranef,
    var = varcomp,
    funvars = parse_functions(reml.out),
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
  tfxd <- stats::terms(fxd)
  
	# Add an intercept manually only if the user requested it
  # *and* there are no other *categorical* fixed effects
  tempmf <- eval(call('model.frame',
                      formula = fxd,
                      data = quote(data)),
                 parent.frame())
  tempt <- stats::terms(tempmf)
  tempc <- attr(tempt, 'dataClasses')[attr(tempt, 'term.labels')]
  any.cat <- any(tempc %in% c('factor', 'ordered'))
	
  if( attr(tfxd, 'intercept') == 1L & !any.cat ) {
	  fxd <- update(fxd, " ~ Intercept + .")
	}
	terms.list$fxd <- attr(stats::terms(fxd), 'term.labels')
	
	
	## Random effects (unstructured)
	rnd <- eval(call$random, parent.frame(2))
  if(!is.null(rnd))
    terms.list$rnd <- attr(stats::terms(rnd), 'term.labels')
	
	## Join fixed and random
	lhs <- deparse(fxd[[2]])
	rhs <- paste(do.call(c, terms.list), collapse = '+')  
	fml <- stats::as.formula(paste(lhs, rhs, sep = '~'), env = parent.frame(2))
	
  # Build Model Frame
  # Use na.pass to allow missing observations which will be handled later
	mfcall <- call('model.frame',
                 formula = fml,
                 data = quote(transform(data, Intercept = 1)),
                 na.action = stats::na.pass)
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


#' Determine a numeric code for missing observations
#' 
#' This function returns a code outside the range of variation of the observed values.
#' 
#' If the range of variation strictly excludes zero, then it is used as the code
#' for missing observations (which is the default code in PROGSF90). Otherwise, 
#' returns 1 - the second power of 10 greater than the maximum observed 
#' magnitude. This ensures a code an order of magnitude greater than the
#' observed values. E.g., for observations in the range -40 -- 28, the code is
#' -999.
#' 
#' @param x a numeric vector.
#' 
#' @examples 
#'   breedR:::pf90_code_missing(rnorm(100))
pf90_code_missing <- function(x) {
  if (min(x, na.rm = TRUE) > 0 || max(x, na.rm = TRUE) < 0) return(0)
  
  ## Second power of 10 greater than maximum observed magnitude
  p10 <- 10^(1+ceiling(log10(max(abs(x), na.rm = TRUE))))
  return(1-p10)
}



#' Default formula for heritability
#' 
#' If all random effects are independent, and there is an additive-genetic
#' effect, computes a default formula in PROGSF90 notation by dividing the
#' genetic variance by the sum of all variance components plus the residual
#' variance.
#' 
#' @return A character vector with one option specification per trait.
#' 
#' @param rglist list of random groups in the parameters of a
#'   \code{\link{progsf90}} object
#' @param traits A character vector with trait names, or NULL for single trait.
#' @param quiet logical. If FALSE, the function issues a message when it fails
#'   to build a formula.
#'   
#' @references 
#'    http://nce.ads.uga.edu/wiki/doku.php?id=readme.aireml#options
pf90_default_heritability <- function (rglist, traits = NULL, quiet = FALSE) {
  
  stopifnot(is.list(rglist))
  
  ## Positions of random effects in the list of random groups
  ranef.idx <- sapply(rglist, function(x) x$pos)
  
  ## trait indices
  tr_idx <- seq_along(traits)
  if (is.null(traits)) tr_idx <- 1
  
  if (length(rglist) > 0 &&                       # there is at least one group
      all(vapply(ranef.idx, length, 1) == 1) &&   # all of them of size 1
      'genetic' %in% names(ranef.idx)) {          # there is a genetic effect
    
    ## Compose a formula term for the variance of random effect x
    ## trait # is 1. vector-friendly
    fterm <- function(x) paste('G', x, x, tr_idx, tr_idx, sep = '_')
    
    ## Additive-genetic variance in the numerator
    ## (Potentially more than one)
    numerator <- fterm(ranef.idx[['genetic']])
    
    ## All variance estimates plus residual variance in the denominator
    trait_component <- cbind(
      matrix(sapply(ranef.idx, fterm), ncol = length(ranef.idx)),
      paste('R', tr_idx, tr_idx, sep = '_')
    )
    denom <- apply(trait_component, 1, paste, collapse = '+')
    
    H2fml <- paste0(numerator, "/(", denom, ")")
    
    H2lbl <- "Heritability"
    if (!is.null(traits)) H2lbl <- paste(H2lbl, traits, sep = ":")
    option.str <- paste('se_covar_function', H2lbl, H2fml)
    
  } else {
    
    option.str <- NULL
    if (!quiet)
      message("Can't compute the heritability formula automatically.")
  }
  
  return(option.str)
}


#' Parse a matrix from a text output robustly
#' 
#' Each row of the matrix is a string. If rows are too long, they can continue
#' in another line. Hence, the number of lines might be a multiple of the number
#' of columns
#' 
#' @param x A character vector with space-separated numbers
#' @param names A character vector with row and column names for the output matrix.
#' @param square logical. Whether to assume that the matrix is square.
#' 
#' If the matrix is not necessarily
parse.txtmat <- function(x, names = NULL, square = TRUE) {
  ## numeric values from the strings, spliting by spaces
  ans <- unname(
    lapply(x, function(x) as.numeric(strsplit(x, ' +')[[1]][-1]))
  )
  
  ## concatenate rows in groups of p
  fix_rows <- function(x, p) {
    grp <- split(x, rep(seq_len(length(x)/p), each = p))
    return(unname(lapply(grp, unlist)))
  }

  row_lengths <- sapply(ans, length)
  nonzero_drl <- diff(row_lengths) != 0L

  ## Handle potential line wrapping
  if( any(nonzero_drl) ) {
    ## row-lengths not constant: e.g. drl = 10 10 2 10 10 2 ...
    ## concatenate lines by period
    period <- head(which(nonzero_drl), 1) + 1
    ans <- fix_rows(ans, period)
    
    ## row lengths should be constant by now
    row_lengths <- sapply(ans, length)
    nonzero_drl <- diff(row_lengths) != 0L
    stopifnot(!any(nonzero_drl))
  }
  
  if ( square && (m <- length(ans)) != (n <- row_lengths[1]) ) {
    ## square matrix not square
    ## check rows are multiple of cols
    if (m %% n != 0) stop("Could not figure out square matrix dimensions.")
    period <- m %/% n
    ans <- fix_rows(ans, period)
  }
  
  ## ensure matrix even if only a number
  ans <- as.matrix(simplify2array(ans))
  
  ## if no sub-names passed, remove all naming
  if( is.null(names) ) names(ans) <- dimnames(ans) <- NULL
  # otherwise, put the given names in both dimensions (covariance matrix)
  else dimnames(ans) <- list(names, names)
  
  return(ans)
}

