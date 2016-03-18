old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

### Test the rendering to pf90 format ###

context("Default heritability")

test_that('pf90_default_heritability() works as expected', {

  ## empty list (no random effects) -> NULL
  rgl <- list()
  expect_null(pf90_default_heritability(rgl))

  ## list with no genetic effect -> NULL
  rgl <- list(a = list(pos = 1), b = list(pos = 2))
  expect_null(pf90_default_heritability(rgl))
  
  ## list with genetic effect, and some group of correlated effects
  ## -> don't know how to compute default herit. -> NULL
  rgl <- list(a = list(pos = 1:2), genetic = list(pos = 3))
  expect_null(pf90_default_heritability(rgl))
  
  ## list with a genetic effect alone -> s_a/(s_a+s_e)
  rgl <- list(genetic = list(pos = 3))
  ans <- "se_covar_function Heritability G_3_3_1_1/(G_3_3_1_1+R_1_1)"
  expect_identical(pf90_default_heritability(rgl), ans)
  
  ## list with a genetic effect and some other effects
  rgl <- list(a = list(pos = 1), genetic = list(pos = 3))
  ans <- "se_covar_function Heritability G_3_3_1_1/(G_1_1_1_1+G_3_3_1_1+R_1_1)"
  expect_identical(pf90_default_heritability(rgl), ans)
})

