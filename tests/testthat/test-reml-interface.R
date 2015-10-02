old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

data(globulus)
ped <- build_pedigree(1:3, data = globulus)
# Test function
fit.model <- function(vi, vigen, random, dat = globulus, ...) {
  res <- 
    try(
      remlf90(fixed   = phe_X ~ gen,
              random  = random,
              var.ini = vi,
              genetic = list(model = 'add_animal', 
                             var.ini = vigen,
                             pedigree = ped,
                             id = 'self'), 
              data    = dat)
      , silent = TRUE)
  res
}

# Test data
testdat <- list(list(vi     = NULL,   # Missing specification: FAIL
                     vigen  = 1,
                     random = ~ bl,
                     expectation = 0),
                list(vi     = list(bl = 1),   # Missing specification of residual: FAIL
                     vigen  = 1,
                     random = ~ bl,
                     expectation = 0),
                list(vi     = list(gg = 1,      # Missing specification of bloc: FAIL
                                   resid = 1),
                     vigen  = 1,
                     random = ~ bl + gg,
                     expectation = 0),
                list(vi     = list(bl = 1,     # Missing specification of gg: FAIL
                                   resid = 1),
                     vigen  = 1,
                     random = ~ bl + gg,
                     expectation = 0),
                list(vi     = list(bl = 1,     # Missing genetic specification: FAIL
                                   resid = 1),
                     vigen  = NULL,
                     random = ~ bl,
                     expectation = 0),
                list(vi     = NULL,              # OK: no specification at all
                     vigen  = NULL,
                     random = ~ bl,
                     expectation = 1),
                list(vi     = list(bl = 1,       # OK
                                   resid = 1),
                     vigen  = 1,
                     random = ~ bl,
                     expectation = 1),
                list(vi     = list(bl = 1,       # OK
                                   gg = 1,
                                   resid = 1),
                     vigen  = 1,
                     random = ~ bl + gg,
                     expectation = 1),
                list(vi     = list(resid = 1),     # OK
                     vigen  = 1,
                     random = NULL,
                     expectation = 1))


#### Context: Variance components specifications ####
context("Variance components specifications")

# reml results
# fit.model(vi=list(resid = 1), vigen=1, random = NULL)
# do.call('fit.model', testdat[[1]])
# do.call('fit.model', testdat[[7]])
reslst <- lapply(testdat, function(x) do.call(fit.model, x))


# Compare expected and true results
run_expectations <- function(m, res) {
  # Check that remlf90 behaves as expected
  test_that("remlf90 requires either full or null variance specifications", {
    ifelse( m$expectation,
            expect_true(!inherits(res, "try-error")),
            expect_true(inherits(res, "try-error")) )
  })
}

for(i in seq_along(testdat)) {
#   cat(i)
  run_expectations(testdat[[i]], reslst[[i]])
}
