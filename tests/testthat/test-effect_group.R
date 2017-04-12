

### Test the effect_group constructor ###
context("effect_group() constructor")


test_that('Valid specs of effect groups pass checks', {
  
  # one-element effect_group
  expect_error(
    effect_group(list(breedr_effect(1)), cov.ini = 1, ntraits = 1),
    NA
  )
  
  expect_error(
    effect_group(list(breedr_effect(1)), cov.ini = diag(2), ntraits = 2),
    NA
  )
  
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
  
  expect_cov_ok <- function(x) {
    eval(bquote(
      expect_error(
        effect_group(eflst, cov.ini = .(x), ntraits = 1),
        NA
      )
    ))
    
    xx <- kronecker(as.matrix(x), diag(2))
    eval(bquote(
      expect_error(
        effect_group(eflst, cov.ini = .(xx), ntraits = 2),
        NA
      )
    ))
  }
  
  for(x in matlst) expect_cov_ok(x)
  
})


test_that('Invalid specs of effect groups are caught by checks', {
  
  x <- list(breedr_effect(1))  # a valid list of (one) effect
  
  # missing x
  # the specific error msg is locale-dependent
  expect_error(effect_group(cov.ini = 1, ntraits = 1))
  
  # missing cov.ini
  # the specific error msg is locale-dependent
  expect_error(effect_group(x, ntraits = 1))
  
  # x not a list
  expect_error(effect_group(1, 1, ntraits = 1), 'is.list.*? is not TRUE')

  # cov.ini not like a matrix
  expect_error(effect_group(x, 'a', ntraits = 1), 
               'is.numeric.*? is not TRUE')

  # element in x not a breedr_effect
  # here x is a list and a breedr_effect, but not its elements.
  expect_error(effect_group(breedr_effect(1), 1, ntraits = 1), 
               'must be of class breedr_effect')
  
  # non coforming dimensions
  expect_error(effect_group(x, diag(1:2), ntraits = 1), 'do not conform')
  expect_error(effect_group(x, cov.ini = 1, ntraits = 2), 'do not conform')
})
