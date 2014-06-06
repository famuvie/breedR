### Test the building of model frames ###

old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

# In particular, check that the intercept attribute is always set to zero,
# and it is manually introduced in the model frame when necessary

context("Model Frame")

# Possible fixed structures, with correct outcome
fxdlst <- list(list(fxd = phe_X ~ 1,
                    int = 1L),
               list(fxd = phe_X ~ sex,
                    int = 1L),
               list(fxd = phe_X ~ 1 + sex,
                    int = 1L),
               list(fxd = phe_X ~ 0 + sex,
                    int = 0L),
               list(fxd = phe_X ~ sex - 1,
                    int = 0L))

# Possible specifications of random effects
rndlst <- list(list(rnd = ~ NULL),
               list(rnd = NULL),
               list(rnd = ~ mum))

# All possible combinations of the fixed and random components
mlst <- unlist(lapply(fxdlst,
               function(x) lapply(rndlst,
                                  function(y) c(x, y))), recursive = FALSE)


# Function that runs checks to each model
run_expectations <- function(m) {
  fc <- call('remlf90',
             fixed  = m$fxd,
             random = m$rnd,
             data   = quote(as.data.frame(m1)))
  
  res <- try(eval(fc), silent = TRUE)

  # The model runs OK
  test_that(paste('The model', fc, 'runs OK'), {
    expect_true(!inherits(res, 'try-error'))
  })
  
  if( !!inherits(res, 'try-error') ){
    mf <- res$mf
    
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
}


# Run expectations
lapply(mlst, run_expectations)
