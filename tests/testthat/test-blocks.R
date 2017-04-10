
#### Context: Blocks infrastructure ####
context("Blocks infrastructure")

test_that("breedr_blocks() constructor gives a list with 5 elements of correct sizes", {
  x.loc <- 1:100
  y.loc <- seq(1000, by = 5, length = 51)
  n.blocks <- 10
  
  full.coord <- expand.grid(x.loc, y.loc)
  coord <- full.coord[sample(nrow(full.coord), nrow(full.coord)/2), ]
  block <- factor(sample(n.blocks, nrow(coord), replace = TRUE))

  result <- breedr_blocks(coord, id = block)
  inc.mat <- model.matrix(result)
  cov.mat <- get_structure(result)

  expect_is(result, c('blocks', 'spatial, random', 'breedr_effect'))
  expect_equal(length(result), 5)
  expect_equal(names(result$param), 'n.blocks')
  expect_equal(n.blocks, result$param$n.blocks)
  expect_equal(nrow(inc.mat), nrow(coord))
  expect_equal(ncol(inc.mat), n.blocks)
  
  # The matrix U should be in sparse format: row col value
  expect_that(cov.mat, is_a('Matrix'))
})
