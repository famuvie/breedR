
data(globulus)
ped <- build_pedigree(1:3, data = globulus)
# Test function
fit.model <- function(vi, vigen, random, dat = globulus, ...) {
  try(
    suppressMessages(
      remlf90(fixed   = phe_X ~ gen,
              random  = random,
              var.ini = vi,
              genetic = list(model = 'add_animal', 
                             var.ini = vigen,
                             pedigree = ped,
                             id = 'self'), 
              data    = dat)
    ),
  silent = TRUE)
}

# Test data
testdat <- list(
  list(
    vi     = NULL,   # Missing specification: FAIL
    vigen  = 1,
    random = ~ bl,
    expectation = 0
  ),
  list(
    vi     = list(bl = 1),   # Missing specification of residual: FAIL
    vigen  = 1,
    random = ~ bl,
    expectation = 0
  ),
  list(
    vi     = list(gg = 1,      # Missing specification of bloc: FAIL
                  resid = 1),
    vigen  = 1,
    random = ~ bl + gg,
    expectation = 0
  ),
  list(
    vi     = list(bl = 1,     # Missing specification of gg: FAIL
                  resid = 1),
    vigen  = 1,
    random = ~ bl + gg,
    expectation = 0
  ),
  list(
    vi     = list(bl = 1,     # Missing genetic specification: FAIL
                  resid = 1),
    vigen  = NULL,
    random = ~ bl,
    expectation = 0
  ),
  list(
    vi     = NULL,              # OK: no specification at all
    vigen  = NULL,
    random = ~ bl,
    expectation = 1
  ),
  list(
    vi     = list(bl = 1,       # OK
                  resid = 1),
    vigen  = 1,
    random = ~ bl,
    expectation = 1
  ),
  list(
    vi     = list(bl = 1,       # OK
                  gg = 1,
                  resid = 1),
    vigen  = 1,
    random = ~ bl + gg,
    expectation = 1
  ),
  list(
    vi     = list(resid = 1),     # OK
    vigen  = 1,
    random = NULL,
    expectation = 1
  )
)


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

#### Context: Multitrait specifications ####
context("Multitrait interface")

## two correlated variables
dim <- 3
Nobs <- 1e4
Nbl <- 50
beta_X <- c(-1, 5, 3)
sample_covar <- function(dim) {
  x <- sample(-2:dim, size = dim**2, replace = TRUE)
  crossprod(matrix(x, nrow = dim))
}
set.seed(123)
S_bl <- sample_covar(dim)   # 5 & 3 & 5 // 22 & 10 // 11
S_resid <- sample_covar(dim)  # 9 & 3 & -3 // 9 & 9 // 14
# diag(1/sqrt(diag(S_bl))) %*% S_bl %*% diag(1/sqrt(diag(S_bl)))

dimnames(S_bl) <- dimnames(S_resid) <- 
  rep(list(paste0("y", seq_len(dim))), 2)

bl_levels <- paste0(
  "bl",
  sprintf(paste0("%0", floor(log10(Nbl)+1), "d"), seq_len(Nbl))
)

testdat <- data.frame(
  X = runif(Nobs),
  breedR.sample.ranef(
    dim, S_bl, Nbl, labels = bl_levels, N = Nobs, vname = 'bl'
  ),
  breedR.sample.ranef(dim, S_resid, Nobs, vname = 'e'))

# var(testdat[, c('e1', 'e2')])  # ~ S_resid

testdat <- 
  transform(testdat,
            y1 = beta_X[1]*X + bl_y1 + e_y1,
            y2 = beta_X[2]*X + bl_y2 + e_y2,
            y3 = beta_X[3]*X + bl_y3 + e_y3)


test_that("Residual variance acurately identified in a fixed-effects model", {
  
  ## All fixed effects (AI fails with 3 traits)
  res <- remlf90(
    cbind(y1, y2, y3) ~ X + bl,
    data = testdat,
    method = "em"
  )
  
  expect_equal(S_resid, res$var$Residual, tol = .01)
})


test_that("Simulated values reasonably recovered using one or more traits", {
  
  ## 1 trait
  res_1 <- remlf90(
    cbind(y1) ~ 0 + X,
    random = ~ bl,
    data = testdat
  )
  
  expect_equal(S_resid[1, 1], res_1$var["Residual", 1], tol = .01)
  expect_equal(S_bl[1, 1], res_1$var["bl", 1], tol = .1)
  expect_equal(beta_X[1], fixef(res_1)$X, tol = .1, check.attributes = FALSE)

  ## 2 trait
  res_2 <- remlf90(
    cbind(y1, y2) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "ai"
  )

  expect_equal(S_resid[-3, -3], res_2$var[["Residual", 1]], tol = .01)
  expect_equal(S_bl[-3, -3], res_2$var[["bl", 1]], tol = 1)
  expect_equal(beta_X[-3], fixef(res_2)$X, tol = .1, check.attributes = FALSE)
  
  ## 3 trait  (with "em" since otherwise does not converge)
  res_3 <- remlf90(
    cbind(y1, y2, y3) ~ 0 + X,
    random = ~ bl,
    data = testdat,
    method = "em"
  )  

  expect_equal(S_resid, res_3$var$Residual, tol = .01)
  expect_equal(S_bl, res_3$var$bl, tol = 1)
  expect_equal(beta_X, fixef(res_3)$X, tol = .01, check.attributes = FALSE)
})



# mf <- model.frame(cbind(V1, V2) ~ 0 + mu, transform(testdat, mu = 1))
# attr(attr(mf, 'terms'), 'term.types') <- list(mu = "fixed")
# eff <- build.effects(mf, genetic = NULL, spatial = NULL, generic = NULL, var.ini = S)
# pf90 <- progsf90(mf, eff, res.var.ini = S)

## Use larix dataset:
## - Two phenotypes: LAS and DOS
## - repeated measurements along 16 years (yr)

inc.mat <- model.matrix(~ 0 + bl, larix)
cov.mat <- diag(nlevels(larix$bl))

test_that("Multitrait model with all kind of effects works as expected", {
  
  fullrun <- function(method, opt = NULL) {
    try(
      suppressMessages(
        remlf90(
          fixed   = cbind(LAS, DOS) ~ rep,
          random  = ~ bl,
          genetic = list(model = 'add_animal',
                         pedigree = larix[, 1:3],
                         id = 'self'),
          spatial = list(model = 'AR',
                         coordinates = larix[, c('x', 'y')],
                         rho = c(.8, .8)),
          generic = list(block = list(inc.mat,
                                      cov.mat)),
          data    = larix,
          method = method,
          progsf90.options = opt
        )
      )
    )
  }
  
  ## make things fast, as I am not looking at numerical results
  res_ai <- fullrun("ai", opt = c("maxrounds 2"))
  ## cannot use it with em as the logfile would not report final estimates
  res_em <- fullrun("em")
  
  fixef_names <- "rep"
  ranef_names <- c("bl", "genetic", "spatial", "block")
  
  
  ## No errors
  expect_false(inherits(res_em, "try-error"))
  expect_false(inherits(res_ai, "try-error"))
  
  ## fixed effect estimates
  expect_identical(names(fixef(res_em)), fixef_names)
  expect_identical(names(fixef(res_ai)), fixef_names)
  
  ## variance component estimates
  expect_identical(names(res_em$var), c(ranef_names, "Residual"))
  expect_identical(names(res_em$var), c(ranef_names, "Residual"))
  
  ## random effect blups
  expect_identical(names(ranef(res_em)), ranef_names)
  expect_identical(names(ranef(res_ai)), ranef_names)
  
})


