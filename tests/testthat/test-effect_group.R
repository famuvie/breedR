old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))


### Test the effect_group constructor ###
context("effect_group() constructor")


test_that('Valid specs of effect groups pass checks', {
  
  # one-element effect_group
  expect_error(effect_group(list(breedr_effect(1)), cov.ini = 1), NA)
  
  # two-element effect_group
  eflst <- list(breedr_effect(1), breedr_effect(2))
  
  matlst <- list(
    diag(1:2),                       # diagonal matrix
    crossprod(matrix(1:4,2)),        # SPD full matrix
    structure(diag(1:2), 
              dimnames = list(1:2, 1:2)),  # named matrix
    structure(diag(1:2), 
              dimnames = list(1:2, 3:4)),  # named matrix with different names
    as.data.frame(diag(1:2)),         # data.frame
    Matrix::Diagonal(2)               # Matrix class
  )
  
  expect_cov_ok <- function(x)
    eval(bquote(expect_error(effect_group(eflst, cov.ini = .(x)), NA)))
  
  for(x in matlst) expect_cov_ok(x)
  
})


test_that('Invalid specs of effect groups are caught by checks', {
  
  x <- list(breedr_effect(1))  # a valid list of (one) effect
  
  # missing x
  expect_error(effect_group(cov.ini = 1))
  
  # missing cov.ini
  expect_error(effect_group(x))
  
  # x not a list
  expect_error(effect_group(1, 1), 'is.list.*? is not TRUE')

  # cov.ini not like a matrix
  expect_error(effect_group(x, 'a'), 
               'is.numeric.*? is not TRUE')

  # element in x not a breedr_effect
  # here x is a list and a breedr_effect, but not its elements.
  expect_error(effect_group(breedr_effect(1), 1), 
               'must be of class breedr_effect')
  
  # cov.ini not symmetric
  expect_error(effect_group(c(x, x), matrix(1:4, 2)), 
               'isSymmetric.*? is not TRUE')
  
  # cov.ini not positive-definite
  expect_error(effect_group(c(x, x), matrix(rep(1, 4), 2)), 
               'all\\(ev > 0\\) is not TRUE')
  
  # non coforming dimensions
  expect_error(effect_group(x, diag(1:2)), 'do not conform')
  
})
