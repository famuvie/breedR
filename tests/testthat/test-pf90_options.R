old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

### Test the interface to PROGSF90 OPTIONS ###

context("PROGSF90 options")

test_that('remlf90 parses progsf90.options correctly', {

  N = 100
  dat <- transform(data.frame(x = runif(N)),
                    y = 1 + 2*x + rnorm(N))
  
  expect_opt <- function(opt, regexp) {
    res <- suppressMessages(
      remlf90(y ~ x,
              data = dat,
              progsf90.options = opt)
    )
    
    for (i in seq_along(opt))
    expect_true(any(grepl(regexp[i], res$reml$output)),
                label = paste('option', opt[i], 'passed correctly'))
  }
  
  ## some additional option
  expect_opt('tol 1d-01', 'tolerance .*? 0\\.1')
  expect_opt('EM-REML 2', 'EM-REML iterations .*? 2')
  
  ## Conflicting option: sol se
  # included
  expect_opt('sol se', 'store solutions and s\\.e\\.')
  # only once
  parameters_file <- readLines(file.path(tempdir(), 'parameters'))
  expect_identical(length(grep('OPTION sol se', parameters_file)), 1L)
  
  ## Conflicting option: missing
  # if explicit, the one set by the user is used
  expect_opt('missing 12345', 'missing observation .*? 12345')
  # only once
  parameters_file <- readLines(file.path(tempdir(), 'parameters'))
  expect_identical(length(grep('OPTION missing', parameters_file)), 1L)

  ## Multiple options, some conflicting, some not
  opts <- c('sol se', 'missing 12345', 'tol 1d-01', 'EM-REML 2')
  expr <- c('store solutions and s\\.e\\.', 
            'missing observation .*? 12345',
            'tolerance .*? 0\\.1',
            'EM-REML iterations .*? 2')
  expect_opt(opts, expr)
  parameters_file <- readLines(file.path(tempdir(), 'parameters'))
  expect_identical(length(grep('OPTION', parameters_file)), length(opts))
})


test_that('remlf90(method = "ai") returns heritability', {
  
  ## Simulate a dataset with a heritability of 2/(1+2+1+1) = 0.4
  dat <- breedR.sample.phenotype(fixed   = c(mu = 10, x = 2),
                                 random = list(u = list(nlevels = 3,
                                                        sigma2  = 1)),
                                 genetic = list(model    = 'add_animal',
                                                Nparents = c(10, 10),
                                                sigma2_a = 2,
                                                check.factorial = FALSE),
                                 spatial = list(model     = 'AR',
                                                grid.size = c(15, 15),
                                                rho       = c(.7, .8),
                                                sigma2_s  = 1),
                                 residual.variance = 1)
  
  res <- remlf90(phenotype ~ 1 + X.x,
                 random = ~ u,
                 genetic = list(model = 'add_animal',
                                pedigree = dat[, 1:3],
                                id = 'self'),
                 spatial = list(model = 'AR',
                                coord = dat[, c('Var1', 'Var2')],
                                rho   = c(.7, .8)),
                 data   = dat)

  expect_true(
    all(
      any(grepl('* SE for function of \\(co\\)variances H2', res$reml$output)),
      any(grepl('H2  - Function: ', res$reml$output))
    )
  )
})
