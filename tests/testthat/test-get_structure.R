### Test the functions for getting the structure matrices ###


#### Context: Extracting structure matrices ####
context("Extracting structure matrices")

## Extracting structure matrices from simple breedR effects

coord <- expand.grid(list(x = seq(1, 100, length = 51),
                          y = seq(1001, 1100, length = 35)),
                     KEEP.OUT.ATTRS = FALSE)

# A splines object
spl <- breedr_splines(coord)
spl.str <- get_structure(spl)

# A generic object (with same structure, but inverted)
inv_spl.str <- solve(spl.str)
gen <- generic(incidence = model.matrix(spl),
               precision = inv_spl.str)
gen.str <- get_structure(gen)

test_that('get_structure() extracts a Matrix', {
  expect_is(spl.str, 'Matrix')
  expect_is(gen.str, 'Matrix')
})


test_that('get_structure() recovers the right structure type', {
  expect_identical(attr(spl.str, 'type'), 'covariance')
  expect_identical(attr(gen.str, 'type'), 'precision')
})


## Extracting structure matrices from groups of effects

eg <- effect_group(list(spl, gen), cov.ini = diag(1,2,2), ntraits = 1)
eg.str <- get_structure(eg)

test_that('get_structure() recovers the common structure in Matrix format', {
  expect_is(eg.str, 'Matrix')
  expect_identical(eg.str, spl.str)
})


## Extracting structure matrices from breedR objects

test_that('get_structure() retrieves an empty list from a model fit without random effects', {

  res <- load_res("fixonly")
  breedr.str <- get_structure(res)
  
  expect_is(breedr.str, 'list')
  expect_equal(breedr.str, list(), check.attributes = FALSE)
})


test_that('get_structure() retrieves a list of structure matrices from a model fit', {

  res <- load_res("ar")
  breedr.str <- get_structure(res)
  
  expect_is(breedr.str, 'list')
  for (i in seq_along(breedr.str)) {
    expect_is(breedr.str[[i]], "Matrix")
  }
  
})
