### Test the functions for getting the structure matrices ###

old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("Extracting structure matrices")

## Extracting structure matrices from simple breedR effects

coord <- expand.grid(list(x = seq(1, 100, length = 51),
                          y = seq(1001, 1100, length = 35)),
                     KEEP.OUT.ATTRS = FALSE)

# A splines object
spl <- splines(coord)
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

eg <- effect_group(list(spl, gen), cov.ini = diag(1,2,2))
eg.str <- get_structure(eg)

test_that('get_structure() recovers the common structure in Matrix format', {
  expect_is(eg.str, 'Matrix')
  expect_identical(eg.str, spl.str)
})