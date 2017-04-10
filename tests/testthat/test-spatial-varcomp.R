
### Scaling of spatial variance components ###
context("Scaling of the spatial variance component")

test_that("The spatial variance component is the characteristic marginal variance of the spatial effect's contribution to the phenotypic variance", {
  
  res <- list(
    blk = load_res("blk"),
    ar = load_res("blk"),
    spl = load_res("blk")
  )
  
  expect_equal(breedR:::gmean(Matrix::diag(vcov(res$blk))), res$blk$var['spatial', 1])
  expect_equal(breedR:::gmean(Matrix::diag(vcov(res$spl))), res$spl$var['spatial', 1])
  expect_equal(breedR:::gmean(Matrix::diag(vcov(res$ar))),  res$ar$var['spatial', 1])
})

