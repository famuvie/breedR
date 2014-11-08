#' Inference with REMLF90
#' 
#' Fits a Linear Mixed Model by Restricted Maximum Likelihood
#' 
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
#' @param var.ini if not \code{NULL}, a named list with one item for each random
#'   component in the model. See 'Details'.
#' @param method either 'ai' or 'em' for Average-Information or 
#'   Expectation-Maximization REML respectively
#' @param breedR.bin character. The local directory where the package binaries 
#'   are stored, or any of 'remote' or 'submit' for remote computing. See 
#'   'Details'.
#' @param debug logical. If \code{TRUE}, the input files for blupf90 programs 
#'   and their output are shown, but results are not parsed.
#'   
#' @details If either \code{genetic} or \code{spatial} are not \code{NULL}, the 
#'   model residuals are assumed to have an additive genetic effects or a 
#'   spatially structured random effect, respectively. In those cases, 
#'   \code{genetic} and \code{spatial} must be lists with named relevant 
#'   parameters.
#'   
#'   \subsection{Genetic effect}{ The available models for the genetic effect 
#'   are \code{add_animal} and \code{competition}. \code{add_animal} stands for 
#'   an additive animal model with a given pedigree. \code{competition} includes
#'   the \emph{direct} additive genetic effect and also a \emph{competition} 
#'   additive genetic effect possibly correlated with the former.
#'   
#'   The minimum elements in the list of the \code{genetic} component are: 
#'   \itemize{\item \code{model} a string, either \code{add_animal} or 
#'   \code{competition} \item \code{pedigree} either an object of class 
#'   \code{pedigree} (see \code{\link{build_pedigree}}) or a \emph{data.frame} 
#'   with exactly three columns identifying the individual, his parent and his 
#'   mother respectively \item \code{id} a string, the name of the variable with
#'   the individual identifyer in the \code{data}}
#'   
#'   Optional common components are:
#'   
#'   \itemize{ \item \code{var.ini} a positive initial value for the variance 
#'   component(s). For a \code{competition} model, the same initial value is 
#'   used for both effects, with a negative initial correlation between them of 
#'   half this value.}
#'   
#'   Finally, for model \code{competition} there are further mandatory and 
#'   optional elements:
#'   
#'   \itemize{ \item mandatory elements \itemize{ \item \code{coord} a matrix, 
#'   list or data.frame with two columns for the rows and columns of the 
#'   observations, respectively. This element is necessary even if duplicated in
#'   a \code{spatial} component. \item \code{competition_decay} a positive 
#'   number. The intensity of competition is weighted by the distance acording 
#'   to \eqn{1/d^\alpha}. This element specifies the exponent \eqn{\alpha} to be
#'   used. Typically \eqn{1} or \eqn{2}. }
#'   
#'   \item optional elements \itemize{ \item \code{pec} Permanent Environmental 
#'   Competition effect. If present, this must be a named list with elements 
#'   \code{present} wich is either \code{TRUE} or \code{FALSE} and (optionally) 
#'   \code{var.ini} specifying the initial variance for this effect. } }
#'   
#'   The Permanent Environmental Competition (\code{pec}) effect is actually 
#'   non-genetic in nature. However, it was included as an option to the 
#'   (genetic) competition effect as it is usually used in conjunction with it. 
#'   }
#'   
#'   
#'   \subsection{Spatial effects}{ The available models for the spatial effect 
#'   are \code{splines} and \code{AR1}. \code{splines} uses a  two-dimensional 
#'   tensor product of B-splines to represent the smooth spatially structured 
#'   effect (Cappa and Cantet, 2007). \code{AR1} uses a kronecker product of 
#'   autoregressive models for the rows and columns.
#'   
#'   In both cases, the minimum necessary components in the list are \itemize{ 
#'   \item \code{model} a string, either \code{splines} or \code{AR} \item 
#'   \code{coord} a matrix, list or data.frame with two columns for the rows and
#'   columns respectively. }
#'   
#'   Optional common components are \itemize{ \item \code{var.ini} a positive 
#'   initial value for the variance component }
#'   
#'   Finally, optional model-dependent components are \itemize{ \item For model 
#'   \code{splines} \itemize{ \item \code{n.knots} a vector of two integers with
#'   the number of \emph{internal} knots for the rows and columns of the spline 
#'   design } \item For model \code{AR} \itemize{ \item \code{rho} a vector of 
#'   two numbers strictly between -1 and 1 with the autoregressive parameters 
#'   for the rows and columns respectively. Alternatively, a matrix or 
#'   data.frame with two columns where every row contain a combination of 
#'   autoregressive parameters to be tried. } }
#'   
#'   The \emph{internal} knots cover the region with observations at regular 
#'   intervals (in each dimension). For the splines design, three additional 
#'   knots are automatically added before the first internal knot, and other 
#'   tree after the last one, in each dimension. As a result, if \code{n.knots =
#'   c(n1, n2)}, then the final number of parameters of the splines model is 
#'   \eqn{(n1 + 2)(n2 + 2)}.
#'   
#'   If \code{n.knots} is omitted, a sensible number of knots for each dimension
#'   is computed based on the number of observations in the experiment. See the 
#'   internal function \code{\link{determine.n.knots}}. This is the default 
#'   function that computes the default number of knots, but you can provide an 
#'   alternative default function through \code{\link{breedR.setOption}}.
#'   
#'   Due to limitations of the REML backend, we can only fit models with 
#'   \emph{fixed} autoregressive parameters. \code{remlf90()} will fit as many 
#'   models as rows in \code{rho}, and return the results of the most likely. It
#'   will also return the list of log-likelihoods for each tried combination of 
#'   autoregressive parameters, in the component \code{rho} of the \code{reml} 
#'   object. This is useful for visualization, and further refinement of the 
#'   search for appropriate parameters.
#'   
#'   If any of the values in either column of \code{rho} is \code{NA}, then a 
#'   default set of values for the corresponding dimension will be set. See 
#'   \code{breedR.getOption('ar.eval')}, for the current defaults. You can set 
#'   your own defaults with \code{\link{breedR.setOption}}. Each will be 
#'   combined with every other value in the other column.
#'   
#'   Omitting the specification of \code{rho} is equivalent to \code{rho = c(NA,
#'   NA)}. }
#'   
#'   \subsection{Intercept}{ An intercept is automatically introduced in the 
#'   model provided the user doesn't explicitly prevents it by using \code{0} or
#'   \code{-1} in the \code{fixed} formula (as conventional in \code{R}), 
#'   \emph{and} there are no other categorical covariates in \code{fixed}. The 
#'   latter condition is actually a limitation of (ai)remlf90 backends, which 
#'   would in any case return an estimate for each level of the categorical 
#'   covariates while returning 0 for the intercept. It does not allow 
#'   alternative parameterizations. }
#'   
#'   
#'   \subsection{Initial variances}{ Initial variance components can also be 
#'   specified through an additional argument \code{var.ini}. You can either use
#'   default initial values for the variance components (see 
#'   \code{?breedR.options}) or specify custom values for \emph{each} and 
#'   \emph{all} variance components in the model. In this case, \code{var.ini} 
#'   must be a named list with one element for each term in \code{random} with 
#'   matching names, plus one last element named \code{residual} for the initial
#'   residual variance. Furthermore if there are \code{genetic} or 
#'   \code{spatial} effects, they must as well include a numeric element 
#'   \code{var.ini} with the initial variance component specification for the 
#'   corresponding effect. }
#'   
#'   
#'   \subsection{Inference method}{ AI-REML is usually faster than EM-REML, and 
#'   it provides more results. Namely, standard errors of the variance 
#'   components estimates, and covariances as well. On the other hand, is less 
#'   robust than EM-REML and it usually gives extreme results when used with the
#'   splines spatial model (as in \code{spatial = list(model ='splines', ...)}).
#'   
#'   Even when an effect accounts for no variance at all, EM-REML will always 
#'   estimate a positive variance which will be determined by the starting 
#'   value. If AI-REML does not converge but EM-REML does with the same dataset 
#'   and model, re-run EM-REML with a small starting value for the effect. If 
#'   the estimate does not change, it is likely that there is no variance. }
#'   
#'   \subsection{Remote computing}{ If \code{breedR.bin = 'remote'}, the REML 
#'   program will be run remotely and the results will be automatically 
#'   transferred back automatically. While if \code{breedR.bin = 'submit'} the 
#'   job will be submitted to the server, and the job-id and other relevant 
#'   information about the model will be returned instantly. The returned object
#'   can be used to retrieve the results or check the status of the job. Several
#'   jobs can be submitted in parallel. See \code{?remote} to learn how to 
#'   configure breedR for remote computing, and how to manage submitted jobs.}
#'   
#' @return An object of class 'remlf90' that can be further questioned by 
#'   \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{fitted}}, etc.
#' @seealso \code{\link[pedigreemm]{pedigree}}
#' @references progsf90 wiki page: \url{http://nce.ads.uga.edu/wiki/doku.php}
#'   
#'   E. P. Cappa and R. J. C. Cantet (2007). Bayesian estimation of a surface to
#'   account for a spatial trend using penalized splines in an individual-tree 
#'   mixed model. \emph{Canadian Journal of Forest Research} 
#'   \strong{37}(12):2677-2688.
#'   
#' @examples
#' ## Linear model
#' n <- 1e3
#' dat <- transform(data.frame(x = runif(n)),
#'                  y = 1 + 2*x + rnorm(n, sd = sqrt(3)))
#' res.lm <- remlf90(fixed = y ~ x, data = dat)
#' summary(res.lm)
#' 
#' ## Animal model
#' ped <- build_pedigree(c('self', 'dad', 'mum'),
#'                       data = as.data.frame(m1))
#' res.am <- remlf90(fixed   = phe_X ~ sex,
#'                   genetic = list(model    = 'add_animal',
#'                                  pedigree = ped,
#'                                  id       = 'self'),
#'                   data    = as.data.frame(m1))
#'                   
#' \dontrun{
#' ## Same model with specification of initial variances
#' res.am <- remlf90(fixed   = phe_X ~ sex,
#'                   genetic = list(model    = 'add_animal',
#'                                  pedigree = ped,
#'                                  id       = 'self',
#'                                  var.ini  = 1),
#'                   data    = as.data.frame(m1),
#'                   var.ini = list(resid = 1))
#'  
#' ## Animal-spatial models
#' gen.globulus <- list(model    = 'add_animal',
#'                      pedigree = globulus[, 1:3],
#'                      id       = 'self')
#' res.bm <- remlf90(fixed   = phe_X ~ gg,
#'                   genetic = gen.globulus,
#'                   spatial = list(model = 'blocks', 
#'                                  coord = globulus[, c('x','y')],
#'                                  id    = 'bl'),
#'                   data    = globulus)
#'                   
#' res.am <- remlf90(fixed   = phe_X ~ gg,
#'                   genetic = gen.globulus,
#'                   spatial = list(model = 'AR', 
#'                                  coord = globulus[, c('x','y')],
#'                                  rho   = c(.85, .8)),
#'                   data    = globulus)
#' 
#' res.sm <- remlf90(fixed   = phe_X ~ gg,
#'                   genetic = gen.globulus,
#'                   spatial = list(model = 'splines', 
#'                                  coord = globulus[, c('x','y')],
#'                                  n.knots = c(5, 5)),
#'                   data    = globulus,
#'                   method  = 'em')   # Necessary for splines models!!!
#' 
#' 
#' ## Competition models
#' 
#' # This may take some minutes...
#' res.cm <- remlf90(fixed   = phe_X ~ 1,
#'                  genetic = list(model = 'competition',
#'                                 pedigree = globulus[, 1:3],
#'                                 id = 'self',
#'                                 coord = globulus[, c('x','y')],
#'                                 competition_decay = 1,
#'                                 pec = list(present = TRUE)),
#'                  data = globulus)
#' }
#' 
#' @export
remlf90 <- function(fixed, 
                    random = NULL,
                    genetic = NULL,
                    spatial = NULL,
                    data, 
                    var.ini = NULL,
                    method = c('ai', 'em'),
                    breedR.bin = breedR.getOption("breedR.bin"),
                    debug = FALSE) {
  
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

  

  ###  Call
  mc <- mcout <- match.call()
  
  ### Checks
  if ( missing(fixed) | missing(data) ) { 
    stop("Usage: remlf90(fixed, data, ...); see ?remlf90\n")
  }
  if ( !inherits(fixed, "formula") ) { 
    stop("'fixed' should be a formula\n") 
  }
  if ( attr(terms(fixed), 'response') != 1L ) {
    stop("There is no response in the 'fixed' argument\n")
  }
  if ( !is.null(random) ) {
    check.random = TRUE
    if ( !inherits(random, "formula") )
      check.random = FALSE
    else if( attr(terms(random), 'response') != 0L )
      check.random = FALSE
    if( !check.random )
      stop("random should be a response-less formula\n")
  }
  # Initial variances specification
  # Either all initial variances specified, or no specification at all
  random.terms <- switch( is.null(random) + 1,
                          c(attr(terms(random), 'term.labels'), 'residuals'),
                          'residuals')
  check.var.ini <- function(eff) {
    ans = FALSE
    if( is.null(eff) ) ans = NA
    else if( !is.null(eff$var.ini) ){
      if( is.numeric(eff$var.ini) & eff$var.ini > 0 ) ans = TRUE
      else stop(paste('var.ini must be > 0 in the',
                      substitute(eff), 'effect.\n'))
    }
    ans
  }
  var.ini.checks <- c(random  = FALSE,
                      genetic = check.var.ini(genetic),
                      pec     = check.var.ini(genetic$pec),
                      spatial = check.var.ini(spatial))
  
  if( !missing(var.ini) ) {
    if( !is.null(var.ini) ) {
      # normalize names
      names(var.ini) <- match.arg(tolower(names(var.ini)),
                                  random.terms,
                                  several.ok = TRUE)
      # check that all the required variances are given
      if( setequal(names(var.ini), random.terms) ) {
        if( all(sapply(var.ini, is.numeric)) & all(var.ini > 0) )
          var.ini.checks[1] = TRUE
        else stop('Initial variances in var.ini must be > 0.\n')
      } else stop('Some initial variances missing in var.ini.\n')
    }
  }
  if( !any(all(var.ini.checks, na.rm = TRUE), all(!var.ini.checks, na.rm = TRUE)) )
    stop('Some initial variances missing. Please check.')
  # In the case of no specification, complete with defaults and Warn
  if( all(!var.ini.checks, na.rm = TRUE) ) {
    warning(paste('No specification of initial variances.\n',
                  '\tUsing default value of',
                  breedR.getOption('default.initial.variance'),
                  'for all variance components.\n',
                  '\tSee ?breedR.getOption.\n'))
    var.ini <- as.list(rep(breedR.getOption('default.initial.variance'),
                           length(random.terms)))
    names(var.ini) <- random.terms
    if( !is.na(var.ini.checks['genetic']) ){
      genetic$var.ini <- breedR.getOption('default.initial.variance')
    }
    if( !is.na(var.ini.checks['spatial']) )
      spatial$var.ini <- breedR.getOption('default.initial.variance')
  }
  
  
  ### Parse arguments
  method <- tolower(method)
  method <- match.arg(method)
  
    

  # Builds model frame by joining the fixed and random terms
  # and translating the intercept (if appropriate) to a fake covariate
  # Add an additional 'term.types' attribute within 'terms'
  # indicating whether the term is 'fixed' or 'random'
	# progsf90 don't allow for custom model parameterizations
	# and they don't use intercepts
  mf <- build.mf(mc)
  mt <- attr(mf, 'terms')

  
  # Genetic effect
  if( !is.null(genetic) ) {
    genetic$model <- match.arg(genetic$model,
                               choices = c('add_animal', 'competition'))

    if( !all(check_pedigree(genetic$pedigree)) ) {
      genetic$pedigree <- build_pedigree(1:3, data = genetic$pedigree)
    }
    
    if( length(genetic$id)==1 ) {
      genetic$id <- data[, genetic$id]
    }
    
    if( genetic$model == 'competition' ) {
      genetic$var.ini <- diag(genetic$var.ini, 2)
      genetic$var.ini[1,2] <- genetic$var.ini[2,1] <- -genetic$var.ini[1,1]/2
      
      # Specification of Permanent Environmental Effect
      if( !exists('pec', genetic) ) {
        # If pec was not listed, put it as not present
        genetic$pec <- list(present = FALSE)
      } else {
        # If it was listed, mark it as present
        if( !exists('present', genetic$pec) ) {
          genetic$pec$present <- TRUE
        }
        # ... and check the initial variance
        if( !exists('var.ini', genetic$pec) ) {
          genetic$pec$var.ini <- breedR.getOption('default.initial.variance')
        }
      }
    }
  }
  
  # Spatial effect
  if( !is.null(spatial) ) {
    spatial$model <- match.arg(spatial$model,
                               choices = c('splines', 'AR', 'blocks'))
    
    # If blocks model, include the values of the relevant covariate
    if( spatial$model == "blocks" ) {
      
      # Check that block id is found in data
      if( !spatial$id %in% names(data) ) {
        stop(paste('Block id variable', spatial$id,
                   'not found in data'))
      }

      spatial$id <- data[, spatial$id]
      
      # Only factors make sense for blocks
      # If it is already a factor, it may have
      # unobserved levels. Otherwise, make it a factor.
      if( !is.factor(spatial$id) )
        spatial$id <- as.factor(spatial$id)
      
    }
    
    # If AR model without rho specified
    # we need to fit it with several fixed rho's
    # and return the most likely
    # TODO: It would be nice if we didn't need to recompute Q each time
    if( spatial$model == 'AR' ) {
      if( is.null(spatial$rho) ) spatial$rho <- matrix(c(NA, NA), 1, 2)
      if( any(is.na(spatial$rho)) ) {
        # Evaluation values for rho
        rho.grid <- build.AR.rho.grid(spatial$rho)
      } else {
        rho.grid <- spatial$rho
      }
      
      if( !is.null(nrow(rho.grid)) ) {
        
        # Results conditional on rho
        eval.rho <- function(rho, mc) {
          mc$spatial$rho <- rho
          ## Bloody function environments and scoping rules!
          mc$spatial$coord <- spatial$coord
          mc$data <- quote(data)
          suppressWarnings(eval(mc))   # Avoid multiple redundant warnings
                                       # about initial variances.
        }
        #         test <- eval.rho(mc, c(.5, .5))
        ans.rho <- apply(rho.grid, 1, eval.rho, mc)
        # Interpolate results
        loglik.rho <- transform(rho.grid,
                                loglik = suppressWarnings(sapply(ans.rho, logLik)))
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
  #   parameter.file.path <- file.path(tmpdir, 'parameters')
  if(!is.null(genetic)) genetic$tempfile <- file.path(tmpdir, 'pedigree')
  if(!is.null(spatial)) spatial$tempfile <- file.path(tmpdir, 'spatial')
  
 
  # Build a list of parameters and information for each effect
  effects <- build.effects(mf, genetic, spatial, var.ini)
  
  # Generate progsf90 parameters
  # TODO: Memory efficiency. At this point there are three copies of the 
  # dataset. One in data, one in mf (only needed variables)
  # and yet one more in pf90. This is a potential problem with large datasets.
  pf90 <- progsf90(mf, effects, opt = c("sol se"), res.var.ini = var.ini$resid)
  
  # Write progsf90 files
  write.progsf90(pf90, dir = tmpdir)

  
  # Change to temporal directory to avoid specification of long paths
  # Avoids Issue #1
  cdir <- setwd(tmpdir)
  on.exit(setwd(cdir))
  
  ## Determine the breedR program to use (either local or remote)
  remote = FALSE
  submit = FALSE
  submit.id = ""
  if ( tolower(breedR.bin) == "remote" || tolower(breedR.bin) == "submit" ) {
    remote = TRUE
    submit.id = paste(gsub("[ :]", "-", date()), "---", as.integer(runif(1,min=1E8,max=1E9-1)), sep="")
    remote.bin = breedR.getOption('remote.bin')
    if( remote.bin == "path_to/breedR/bin/linux" ) {
      stop('breedR is not configured for remote computing. See ?breedR.options')
    }
    
    if( breedR.os('windows') ) {
      if( !breedR.cygwin.check() ) {
        stop(paste("Cannot find the CYGWIN installation:", breedR.getOption("cygwin")))
      }
      
      # Make sure the binaries are accesible
      breedR.cygwin.setPATH()
    }
    
    breedR.call = switch(method,
                         ai = file.path(remote.bin, 'airemlf90'),
                         em = file.path(remote.bin, 'remlf90'))
    
    # Run either breedR.remote or breedR.submit
    if ( tolower(breedR.bin) == "remote" ) {
      ldir <- breedR.remote(submit.id, breedR.call)

      if( !identical(normalizePath(ldir), normalizePath(tmpdir)) ) stop('This should not happen')
      reml.out <- readLines(file.path(ldir, 'LOG'))
    }
    
    if ( tolower(breedR.bin) == "submit" ) {
      submit = TRUE
      reml.out <- breedR.submit(submit.id, breedR.call)
    }
  } else {
    breedR.call = switch(method,
                         ai = file.path(breedR.bin, 'airemlf90'),
                         em = file.path(breedR.bin, 'remlf90'))

    reml.out <- system2(breedR.call, 
                        input  = 'parameters',
                        stdout = ifelse(debug, '', TRUE))
  }
  
  
  
  if( !debug ) {
    # Error catching
    stopifnot(is.null(attr(reml.out, 'status')))
    
    if( !submit ) {
      # Parse solutions
      ans <- parse_results(file.path(tmpdir, 'solutions'), effects, mf, reml.out, method, mcout)
    } else {
      # Submitted job (solutions are parsed later with breedR.qget)
      ans <- list(id = submit.id,
                  effects = effects,
                  mf = mf,
                  method = method,
                  mcout = mcout)
    }
    
    class(ans) <- c('breedR', 'remlf90')  # Update to merMod in newest version of lme4 (?)
  } else {
    file.show('parameters')
    ans = NULL
  }
  
  return(ans)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Internal methods  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#



#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Interface methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @describeIn get_pedigree Get the pedigree from a remlf90 object
#' @export
get_pedigree.remlf90 <- function(x, ...) {
  ped <- x$effects$genetic$ped
  if( !is.null(ped) ) {
    map <- attr(ped, 'map')
    ped <- with(ped,
                pedigreemm::pedigree(sire=sire, dam=dam, label=self))
    attr(ped, 'map') <- map
  }
  return(ped)
}

#' @export
coef.remlf90 <- function(object, ...) { 
  c(fixef(object), ranef(object))
}

#' @export
extractAIC.remlf90 <- function(fit, scale, k, ...) {
  return(fit$fit$AIC)
}
  
#' @method fitted remlf90
#' @export
fitted.remlf90 <- function (object, ...) {
  object$mu
}
  
#' @importFrom nlme fixef
#' @export fixef
#' @export
fixef.remlf90 <- function (object, ...) {
      object$fixed
}

#' @method logLik remlf90
#' @export
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

#' @method model.frame remlf90
#' @export
model.frame.remlf90 <- function (formula, ...) {
  formula$mf
}

#' @method model.matrix remlf90
#' @export
model.matrix.remlf90 <- function (object, ...) {
  object$mm
}

#' @importFrom stats nobs
#' @method nobs remlf90
#' @export
nobs.remlf90 <- function (object, ...) {
  nrow(as.matrix(object$y))
}


# sp::coordinates() is an S4 function
# Register the S3 class 'breedR' as an S4 class
setOldClass('breedR')
setMethod('coordinates', signature = 'breedR', 
          function(obj, ...) {
            if( !obj$components$spatial ) {
              if( !obj$components$pedigree ) {
                stop("This breedR object has no spatial structure.\n")
              } else {
                stopifnot(exists('coord', obj$effects$genetic$gen))
                return(obj$effects$genetic$gen$coord)
              }
            }
            return(obj$effects$spatial$sp$coord)
          }
)

setMethod('coordinates<-', signature = 'breedR', 
          function(object, value) {
            # sp::coordinates() performs some sanity checks
            # like coords to be positive integers, not NA, etc.
            cc <- as.data.frame(sp::coordinates(value))
            if( !object$components$spatial ) {
              object$effects$spatial <- list(name = 'none',
                                             sp = list(coord = cc))
              # Now it is a "spatial" object, although there is no
              # spatial model
              object$components$spatial = TRUE
            } else {
              object$effects$spatial$sp$coord = cc
            }
            return(object)
          }
)


#' Spatial plot of a model's fit components
#' 
#' Plots the predicted values of the component effects of the phenotype.
#' 
#' @param x A \code{breedR} object
#' @param type Character. Which component is to be represented in the map.
#' @param z Optional. A numeric vector to be plotted with respect to the spatial
#'   coordinates. Overrides \code{type}.
#' @param ... Further layers passed to \code{\link[ggplot2]{ggplot}}.
#'   
#' @method plot remlf90
#' @export
plot.remlf90 <- function (x, type = c('phenotype', 'fitted', 'spatial', 'fullspatial', 'residuals'), z = NULL, ...) {
  
  type = match.arg(type)
  
  coord <- try(coordinates(x), silent = TRUE)
  if( inherits(coord, "try-error") ) {
    stop(paste('Missing spatial structure. Use coordinates(',
               deparse(substitute(x)),
               ') <- coord', sep = ''))
  }
  
  # Argument z is used for plotting a custom spatial variable
  if( !is.null(z) ) {
    z <- as.vector(z)
    if( !is.numeric(z) | length(z) != nrow(coord) )
      stop(paste("'z' must be a vector of length", nrow(coord), ".\n"))
    
    # Determine the type of scale, depending on whether the vector includes 0
    if( min(z) < 0 & max(z) > 0 ) sc = 'div'
    else sc = 'seq'
    
    p <- spatial.plot(data.frame(coord, z = z), scale = sc)
    
  } else {
    
    if(type == 'phenotype') {
      
      mf <- x$mf
      mt <- attr(mf, 'terms')
      resp.idx <- attr(mt, 'response')
      if( resp.idx == 0L ) stop('There is no response in the model frame')
      resp <- with(mf, eval(attr(mt, 'variables'))[[resp.idx]])
      
      p <- spatial.plot(data.frame(coord, z = resp), scale = 'seq')
    }
    
    
    if(type == 'fitted') {
      
      p <- spatial.plot(data.frame(coord, z = fitted(x)), scale = 'seq')
    }
    
    
    if(type == 'spatial' | type == 'fullspatial'){
      if( x$components$spatial ) {
        
        if( type == 'spatial' ) {
          spdat <- data.frame(coord,
                              z = x$spatial$fit$z,
                              model = x$call$spatial$model)
        } else {
          spdat <- data.frame(x$spatial$prediction,
                              model = x$call$spatial$model)
        }
        
        p <- spatial.plot(spdat, scale = 'div') + facet_wrap(~ model)
      } else stop('This model has no spatial effect')
    }
    
    if(type == 'residuals') {
      p <- spatial.plot(data.frame(coord, z = residuals(x)), scale = 'div')
    }
  }
  
  # Further args
  if( !missing(...) ) {
    p <- p + ...
  }
  
  p
}

# predict.remlf90 <- function (object, ...) {
#   
# }
# 
# print.remlf90 <- function (object, ...) {
#   
# }

#' @importFrom nlme ranef
#' @export ranef
#' @export
ranef.remlf90 <- function (object, ...) {
  object$ranef
}


#' Covariance matrix of a fitted remlf90 object
#' 
#' Returns the variance-covariance matrix of the specified random effect.
#' For the moment, only the spatial effect is of interest.
#' 
#' @param object a fitted model of class \code{remlf90}
#' @param effect the structured random effect of interest
#' @param ... Not used.
#' 
#' @method vcov remlf90
#' @export
vcov.remlf90 <- function (object, effect = 'spatial', ...) {
  
  effect <- match.arg(effect)
  if( !exists(effect, envir = as.environment(object$effects)) )
    stop(paste('There is no', effect, 'effect in this object.\n'))
  
  # Underlying covariance matrix
  # (only the lower triangle)
  Usp <- object$effects[[effect]]$sp$U
  dimU <- c(max(Usp[, 1]), max(Usp[, 2]))
  stopifnot(identical(dimU[1], dimU[2]))
  
  U <- Matrix::sparseMatrix(i = Usp[, 1], j = Usp[, 2], x = Usp[, 3], symmetric = TRUE)
  
  # need to invert?
  if(object$effects[[effect]]$sp$Utype == 'precision')
    U <- solve(U)
  
  # Incidence matrix
  # It can be a vector if it is an identity reordered (e.g. block or AR cases)
  Bsp <- as.matrix(object$effects[[effect]]$sp$B)
  
  if( ncol(Bsp) == 1 ) {
    dimB = nrow(Bsp)
    B <- Matrix::spMatrix(dimB, dimU[1], i = 1:dimB, j = Bsp[, 1], x = rep(1, dimB))
  } else {
    B <- Matrix::Matrix(Bsp, sparse = TRUE)
  }

  # Scaling parameter
  sigma2 <- object$var[effect, 1]
  
  V <- suppressMessages(sigma2 * B %*% U %*% Matrix::t(B))
  return(V)
}


#' @method residuals remlf90
#' @export
residuals.remlf90 <- function (object, ...) {
  # TODO: to be used when na.action is included in remlf90
  #   naresid(object$na.action, res)
  # TODO: add a scale parameter to allow studentization
  #     first need to determine the right sigma for each residual
  object$residuals
}

#' @method summary remlf90
#' @export
summary.remlf90 <- function(object, ...) {
  
  # If this is a submitted job, return the corresponding qstat object
  # instead of a summary.remlf90 object
  if( exists('id', object) ) {
    # Submitted Job
    return(breedR.qstat(object))
  } 
  
  ans <- object
  
  # Literal description of the model
  effects <- paste(names(ans$components), sep=' and ')
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
#' @method print remlf90
#' @export
print.remlf90 <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  # If this is a submitted job, return the corresponding qstat object
  # instead of a summary.remlf90 object
  if( exists('id', x) ) {
    # Submitted Job
    return(print(breedR.qstat(x)))
  } 
  
  cat(x$model.description, '\n')
  if(!is.null(x$call$formula))
    cat("Formula:", x$formula,"\n")
  if(!is.null(x$call$data))
    cat("   Data:", deparse(x$call$data),"\n")
  if(!is.null(x$call$subset))
    cat(" Subset:", x$call$subset,"\n")
  print(x$model.fit, digits = digits)
  
  if( x$components$spatial & !is.null(x$spatial$name)) {
    switch(x$spatial$name,
           AR = cat(paste("\nAutoregressive parameters for rows and columns: (",
                          paste(x$spatial$model$param, collapse = ', '),
                          ")\n", sep = '')),
           splines = cat(paste("\nNumber of inner knots for rows and columns: (",
                               paste(x$spatial$model$param, collapse =', '),
                               ")\n", sep = ''))
    )
  }
  
  invisible(x)
}



## This is modeled a bit after  print.summary.lm :
#' @method print summary.remlf90
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
  
  if( x$components$spatial & !is.null(x$spatial$name)) {
    switch(x$spatial$name,
           AR = cat(paste("\nAutoregressive parameters for rows and columns: (",
                    paste(x$spatial$model$param, collapse = ', '),
                    ")\n", sep = '')),
           splines = cat(paste("\nNumber of inner knots for rows and columns: (",
                         paste(x$spatial$model$param, collapse =', '),
                         ")\n", sep = ''))
    )
  }
  
  cat("\nVariance components:\n")
  print(x$var, quote = FALSE, digits = digits, ...)
  
  cat('\nFixed effects:\n')
  printCoefmat(x$coefficients)
  invisible(x)
}
