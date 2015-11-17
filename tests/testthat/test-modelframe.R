### Test the building of model frames ###

old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

# In particular, check that the intercept attribute is always set to zero,
# and it is manually introduced in the model frame when necessary

context("Model Frame")

data <- transform(as.data.frame(m1),
                  mum2 = mum,
                  gen = factor(sample(1:2, size = 1600, replace = TRUE)))

# Possible fixed structures, with correct outcome
fxdlst <- list(list(fxd = phe_X ~ 1,
                    int = 1L),
               list(fxd = phe_X ~ sex,
                    int = 0L),
               list(fxd = phe_X ~ 1 + sex + gen,
                    int = 0L),
               list(fxd = phe_X ~ 0 + sex,
                    int = 0L),
               list(fxd = phe_X ~ sex + gen - 1,
                    int = 0L))

# Possible specifications of random effects
rndlst <- list(list(rnd = NULL),
               list(rnd = ~ mum),
               list(rnd = ~ mum + mum2))

# All possible combinations of the fixed and random components
mlst <- unlist(lapply(fxdlst,
               function(x) lapply(rndlst,
                                  function(y) c(x, y))), recursive = FALSE)


# Function that runs checks to each model spec
run_expectations <- function(m) {
  fc <- call('remlf90',
             fixed  = m$fxd,
             random = m$rnd,
             data   = quote(as.data.frame(m1)))
  
  mf <- try(build.mf(fc))

  ## The model frame builds OK
  test_that(paste('The model', deparse(fc), 'runs OK'), {
    expect_true(!inherits(mf, 'try-error'))
  })
  
  ## Handling of intercept
  if( !inherits(mf, 'try-error') ){
    
    # The intercept attr. should always be 0 to avoid problems later
    # with the model matrix. The intercept is manually introduced.
    test_that(paste('The intercept attribute is set to 0 for model', fc), {
      expect_identical(attr(attr(mf, 'terms'), 'intercept'), 0L)
    })
    
    # If there was an intercept either implicit or explicit in the formula
    # There should be an "Intercept" covariate full of 1s, and only in that case.
    test_that(paste('An intercept is manually added to model', fc, 'if and only if it is required'), {
      expect_true(m$int == 1 & exists('Intercept', mf) | m$int == 0 & !exists('Intercept', mf))
    })
  }
  
  ## Building the effects list
  if( !inherits(mf, 'try-error') ){
    
    fixef <- attr(terms(fc$fixed), 'term.labels')    
    if (m$int == 1) fixef <- c('Intercept', fixef)
    
    if(is.null(fc$random)) {
      ranef <- vector('character')
    } else ranef <- attr(terms(fc$random), 'term.labels')

    eff <- build.effects(mf = mf,
                         genetic = NULL,
                         spatial = NULL,
                         generic = NULL,
                         var.ini = sapply(c(ranef, 'residuals'), function(x) 1))

    
    ## The effects list builds OK
    test_that(paste('build.effects() builds the model', deparse(fc)), {
      expect_true(!inherits(eff, 'try-error'))
    })
    
    if( !inherits(eff, 'try-error') ){
      ## Parse correctly the variables in model frame
      test_that(paste('All variables have been accounted for, in model', fc), {
        expect_identical(names(mf)[-1], names(eff))
      })
      
      test_that(paste('Fixed effects recognised in model', fc), {
        expect_identical(fixef,
                         names(which(sapply(eff, inherits, 'fixed'))))
      })
      
      test_that(paste('Random effects recognised in model', fc), {
        expect_identical(ranef,
                         names(which(sapply(eff, inherits, 'effect_group'))))
      })
    }
  }
    
}


# Run expectations
for (m in mlst) run_expectations(m)


## TODO: 
## - introduce minimal specifications of all special components
##   i.e. add_animal, competition (with and without pec), blocks, ar, splines, generic
## - build all possible combinations with mlst
## - adapt tests

m <- mlst[[15]]

gen_spec <- check_genetic(model = 'competition',
                          pedigree = data[, c('self', 'dad', 'mum')],   # pedigree
                          id = data$self,          # vector
                          coordinates = data[, c('irow', 'icol')],  # matrix
                          pec = TRUE)

sp_spec <- check_spatial(model = 'splines',
                         coordinates = data[, c('irow', 'icol')])

x1 <- list(inc = matrix(1,1600,3), cov = diag(3), var.ini = 6)
x2 <- list(inc = matrix(1:8,1600,2), pre = diag(2), var.ini = 4)
x <- list (a = x1, b = x2)
grc_spec <- check_generic(x)

fc <- call('remlf90',
           fixed  = m$fxd,
           random = m$rnd,
           genetic = gen_spec,
           spatial = sp_spec,
           generic = grc_spec,
           data   = quote(as.data.frame(m1)))

mf <- try(build.mf(fc))


if(is.null(fc$random)) {
  ranef <- vector('character')
} else ranef <- attr(terms(fc$random), 'term.labels')

eff <- try(
  build.effects(mf = mf,
                genetic = fc$genetic,
                spatial = fc$spatial,
                generic = fc$generic,
                var.ini = sapply(c(ranef, 'residuals'), function(x) 1))
)

test_that('Build a full model frame with all componenents', {
  expect_false(inherits(eff, 'try-error'))
})


## Compile progsf90 object
context('Compile progsf90 object')

test_that('Compile a full model frame with all componenents', {
  expect_error(
    progsf90(mf, eff, opt = c("sol se"), res.var.ini = 10), 
    NA
  )
})

test_that('PROGSF90 codes for missing values', {
  ## strictly positive or negative observations
  expect_identical(pf90_code_missing(c(.1, 9, NA)), 0)
  expect_identical(pf90_code_missing(c(-.1, -9, NA)), 0)
  
  ## range include 0 (even in the boundary)
  ## code always an order of magnitude above observations
  expect_identical(pf90_code_missing(c(0, 0.9, NA)), -9)
  expect_identical(pf90_code_missing(c(-.1, 9, NA)), -99)
  expect_identical(pf90_code_missing(c(-11, .1, NA)), -999)
  expect_identical(pf90_code_missing(c(-110, 345, NA)), -9999)
})


test_that('If phenotype excludes 0, use default missing code', {
  ## make sure we have posiive observations
  mf$phe_X <- 1 - min(mf$phe_X) + mf$phe_X
  mf$phe_X[1] <- NA
  res.try <- expect_error(
    pf90 <- progsf90(mf, eff, opt = c("sol se"), res.var.ini = 10),
    NA
  )
  
  if (res.try$passed) {
    ## no explicit 'missing' (using default of 0)
    expect_false(any(grepl('missing', pf90$parameter$options)))
	
	## observation value set at 0
    expect_equal(pf90$data[1, 'phenotype'], 0, check.attributes = FALSE)
  }
})

test_that('If phenotype includes 0, use alternative missing code', {
  mf$phe_X[1] <- NA
  res.try <- 
    expect_error(
      pf90 <- progsf90(mf, eff, opt = c("sol se"), res.var.ini = 10), 
      NA
    )
  
  if (res.try$passed) {
    ## explicit 'missing' 
    missing_code <- pf90_code_missing(mf$phe_X)
    expect_true(paste('missing', missing_code) %in% pf90$parameter$options)
    
    ## observation value set at corresponding missing code
    expect_equal(pf90$data[1, 'phenotype'], 
                 missing_code,
                 check.attributes = FALSE)
  }
})


