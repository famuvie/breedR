### Test the building of model frames ###


# In particular, check that the intercept attribute is always set to zero,
# and it is manually introduced in the model frame when necessary

context("Model Frame infrastructure")

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
  fc <- call('remlf90',   # call only, no eval
             fixed  = m$fxd,
             random = m$rnd,
             data   = quote(as.data.frame(m1)))
  
  mf <- try(build.mf(fc))

  fc_txt <- paste(deparse(fc), collapse = " ")
  
  ## The model frame builds OK
  test_that(paste('The model', fc_txt, 'runs OK'), {
    expect_true(!inherits(mf, 'try-error'))
  })
  
  ## Handling of intercept
  if( !inherits(mf, 'try-error') ){
    
    # The intercept attr. should always be 0 to avoid problems later
    # with the model matrix. The intercept is manually introduced.
    test_that(paste('The intercept attribute is set to 0 for model', fc_txt), {
      expect_identical(attr(attr(mf, 'terms'), 'intercept'), 0L)
    })
    
    # If there was an intercept either implicit or explicit in the formula
    # There should be an "Intercept" covariate full of 1s, and only in that case.
    test_that(paste('An intercept is manually added to model', fc_txt, 'if and only if it is required'), {
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
    test_that(paste('build.effects() builds the model', fc_txt), {
      expect_true(!inherits(eff, 'try-error'))
    })
    
    if( !inherits(eff, 'try-error') ){
      ## Parse correctly the variables in model frame
      test_that(paste('All variables have been accounted for, in model', fc_txt), {
        expect_identical(names(mf)[-1], names(eff))
      })
      
      test_that(paste('Fixed effects recognised in model', fc_txt), {
        expect_identical(fixef,
                         names(which(sapply(eff, inherits, 'fixed'))))
      })
      
      test_that(paste('Random effects recognised in model', fc_txt), {
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
                          pec = TRUE,
                          response = data$phe_X)

sp_spec <- check_spatial(model = 'splines',
                         coordinates = data[, c('irow', 'icol')],
                         response = data$phe_X)

x1 <- list(inc = matrix(1,1600,3), cov = diag(3), var.ini = 6)
x2 <- list(inc = matrix(1:8,1600,2), pre = diag(2), var.ini = 4)
x <- list (a = x1, b = x2)
grc_spec <- check_generic(x, response = data$phe_X)

fc <- call('remlf90',   # call only, no eval
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
    progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10), 
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
  
  expect_error(
    pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10),
    NA
  )
  
  ## no explicit 'missing' (using default of 0)
  expect_false(any(grepl('missing', pf90$parameter$options)))
  
	## observation value set at 0
  expect_equal(pf90$data[1, 'phenotype'], 0, check.attributes = FALSE)

  })

test_that('If phenotype includes 0, use alternative missing code', {
  mf$phe_X[1] <- NA

  expect_error(
    pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10), 
    NA
  )
  
  ## explicit 'missing' 
  missing_code <- pf90_code_missing(mf$phe_X)
  expect_true(paste('missing', missing_code) %in% pf90$parameter$options)
  
  ## observation value set at corresponding missing code
  expect_equal(pf90$data[1, 'phenotype'], 
               missing_code,
               check.attributes = FALSE)
})


test_that('progsf90() forbids missing values in fixed effects', {
  mf$sex[1] <- NA;
  eff <- try(
    build.effects(mf = mf,
                  genetic = fc$genetic,
                  spatial = fc$spatial,
                  generic = fc$generic,
                  var.ini = sapply(c(ranef, 'residuals'), function(x) 1))
  )
  
  expect_error(
    pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10),
    'Missing values in covariates are not allowed'
  )
})


test_that('progsf90() admits missing values in random effects', {
  mf$mum[1] <- NA;
  eff <- build.effects(mf = mf,
                       genetic = fc$genetic,
                       spatial = fc$spatial,
                       generic = fc$generic,
                       var.ini = sapply(c(ranef, 'residuals'), function(x) 1))
  pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10)
  
  # the incidence matrix for the first individual is zero
  expect_true(identical(sum(model.matrix(eff$mum$effects[[1]])[1,]), 0))
  cols <- as.numeric(head(strsplit(pf90$parameter$effects$mum, ' ')[[1]], 1))
  expect_true(identical(unname(pf90$data[1, cols]), rep(0, length(cols))))

  # blocks case
  fc$spatial$model <- 'blocks'
  fc$spatial$id <- mf$mum
  eff <- build.effects(mf = mf,
                       genetic = fc$genetic,
                       spatial = fc$spatial,
                       generic = fc$generic,
                       var.ini = sapply(c(ranef, 'residuals'), function(x) 1))
  pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10)

  # the incidence matrix for the first individual is zero
  expect_true(identical(sum(model.matrix(eff$spatial$effects[[1]])[1,]), 0))
  cols <- as.numeric(head(strsplit(pf90$parameter$effects$spatial, ' ')[[1]], 1))
  expect_true(identical(unname(pf90$data[1, cols]), rep(0, length(cols))))
  
})

test_that('progsf90() admits missing values in coordinates', {
  fc$spatial$coordinates$irow[1] <- NA;
  
  # splines
  eff <- build.effects(mf = mf,
                       genetic = fc$genetic,
                       spatial = fc$spatial,
                       generic = fc$generic,
                       var.ini = sapply(c(ranef, 'residuals'), function(x) 1))
  pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10)
  
  # the incidence matrix for the first individual is zero
  expect_true(identical(sum(model.matrix(eff$spatial$effects$splines)[1,]), 0))
  cols <- as.numeric(sapply(strsplit(pf90$parameter$effects$spatial, ' '),
                            head, 1))
  expect_true(identical(unname(pf90$data[1, cols]), rep(0, length(cols))))
  
  # AR
  fc$spatial$model <- 'AR'
  fc$spatial$rho <- c(.9, .9)
  
  eff <- build.effects(mf = mf,
                       genetic = fc$genetic,
                       spatial = fc$spatial,
                       generic = fc$generic,
                       var.ini = sapply(c(ranef, 'residuals'), function(x) 1))
  pf90 <- progsf90(mf, w = NULL, eff, opt = c("sol se"), res.var.ini = 10)
  
  # the incidence matrix for the first individual is zero
  expect_true(identical(sum(model.matrix(eff$spatial$effects$ar)[1,]), 0))
  
  cols <- as.numeric(head(strsplit(pf90$parameter$effects$spatial, ' ')[[1]], 1))
  expect_true(identical(unname(pf90$data[1, cols]), rep(0, length(cols))))
})


# ## How do lme4 and INLA deal with missing values in effects?
# require(lme4)
# require(INLA)
# mf$sex[1] <- mf$mum[2] <- NA;
# 
# res.lme4 <- lmer(phe_X~ sex + (1|mum), data = mf)
# nrow(model.frame(res.lme4))  # 1598: removes observations
# 
# res.inla <- inla(
#   phe_X~ sex + f(mum, model = 'iid'), 
#   data = mf,
#   # INLA gives an error, unless expand.factor.strategy = 'inla'
#   control.fixed = list(expand.factor.strategy='inla'),
#   # return results for the linear predictor
#   control.predictor = list(compute = TRUE)
# )
# 
# head(res.inla$model.matrix)  # uses 0 in the incidence matrix
# str(res.inla$size.linear.predictor) # does not remove anything
