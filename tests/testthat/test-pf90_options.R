old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

### Test the interface to PROGSF90 OPTIONS ###

context("PROGSF90 options")

test_that('remlf90 passes progsf90.options correctly', {

  N = 100
  dat <- transform(data.frame(x = runif(N)),
                    y = 1 + 2*x + rnorm(N))
  mf <- build.mf()
  
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

