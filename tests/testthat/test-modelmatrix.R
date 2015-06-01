### Test the computation of model matrices ###

old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("Model Matrix")

fullcoord <- expand.grid(list(x = seq(1, 100, length = 51),
                          y = seq(1001, 1100, length = 35)),
                     KEEP.OUT.ATTRS = FALSE)
fullN <- nrow(fullcoord)
rm.idx <- sample(fullN, fullN/10)   
coord <- fullcoord[-rm.idx,]

## TODO

## Test also the auxiliar function matrix.short16
## used to convert format from 16 cols to sparse.


## Incidence matrix for splines objects

# Remove 10% of the grid and build a splines model from it
spl <- breedr_splines(coord)

test_that('model.matrix() works as expected with spline objects', {
  expect_equal(nrow(model.matrix(spl)), nrow(coord))
})

test_that('model.matrix(·, fullgrid=TRUE) works as expected with spline objects', {
  mmsplfg <- model.matrix(spl, fullgrid = TRUE)
  expect_equal(nrow(mmsplfg), fullN)
  expect_equal(nrow(attr(mmsplfg, 'coordinates')), fullN)
})


## Incidence matrix for generic objects
## As they are not 'spatial', they don't have coordinates
## and model.matrix() should treat them as mere breedr_effects

# A generic object
gen.str <- get_structure(spl)
gen.inc <- as(sample(ncol(gen.str), nrow(coord), replace = TRUE), 'indMatrix')
gen <- generic(gen.inc, gen.str)


test_that('model.matrix() works as expected with generic objects', {
  mmgen <- model.matrix(gen)
  expect_equal(nrow(mmgen), nrow(coord))
  expect_identical(mmgen, gen.inc)
})

test_that('model.matrix(·, fullgrid=TRUE) works as expected with generic objects', {
  expect_identical(model.matrix(gen, fullgrid = TRUE), model.matrix(gen))
})


## Extracting incidence matrices from groups of effects
## with identical covariance structure

eg <- effect_group(list(spl, gen), cov.ini = diag(1,2,2))
eg.mm <- model.matrix(eg)

test_that('model.matrix() recovers the set of incidence matrices binded by columns', {
  expect_is(eg.mm, 'Matrix')
  expect_identical(eg.mm, Matrix::cBind(model.matrix(spl), gen.inc))
})
